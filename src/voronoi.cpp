/**
 * @file voronoi.cpp
 * @brief Implementation of Fortune's Algorithm for Voronoi Diagrams
 *
 * This file contains the complete implementation of Fortune's sweep line
 * algorithm. The algorithm works by maintaining a "beach line" of parabolas
 * and processing events as the sweep line moves from top to bottom.
 *
 * Time Complexity: O(n log n)
 * Space Complexity: O(n)
 *
 * @author Jack Woods, Ashleigh Kirkpatrick
 * @date 2025
 */

#include "voronoi.h"

namespace Voronoi
{

    // ============================================================================
    // CONSTRUCTOR AND DESTRUCTOR
    // ============================================================================

    /**
     * @brief Default constructor
     *
     * Initializes the Voronoi diagram with default bounding box and empty
     * data structures.
     */
    VoronoiDiagram::VoronoiDiagram()
        : beachLineRoot(nullptr),
          visibleMinX(0), visibleMaxX(800), visibleMinY(0), visibleMaxY(600),
          boundMinX(-COMPUTE_MARGIN), boundMaxX(800 + COMPUTE_MARGIN),
          boundMinY(-COMPUTE_MARGIN), boundMaxY(600 + COMPUTE_MARGIN),
          sweepLineY(0)
    {
        // Data structures initialized through member initialization list
    }

    /**
     * @brief Destructor
     *
     * Ensures all dynamically allocated memory is properly freed.
     */
    VoronoiDiagram::~VoronoiDiagram()
    {
        cleanup();
    }

    // ============================================================================
    // PUBLIC METHODS - SITE MANAGEMENT
    // ============================================================================

    /**
     * @brief Add a site to the diagram using coordinates
     *
     * Each site will become a region in the final Voronoi diagram.
     * Sites should be added before calling compute().
     *
     * @param x X coordinate of the site
     * @param y Y coordinate of the site
     */
    void VoronoiDiagram::addSite(double x, double y)
    {
        Point p(x, y, static_cast<int>(sites.size()));
        sites.push_back(p);
    }

    /**
     * @brief Add a site using a Point object
     *
     * @param p The point to add (id will be assigned automatically)
     */
    void VoronoiDiagram::addSite(const Point& p)
    {
        Point newPoint = p;
        newPoint.id = static_cast<int>(sites.size());
        sites.push_back(newPoint);
    }

    /**
     * @brief Clear all data and reset to initial state
     *
     * Removes all sites, edges, and computed data. Must be called
     * if you want to reuse the object for a new diagram.
     */
    void VoronoiDiagram::clear()
    {
        cleanup();
        sites.clear();
        edges.clear();
        beachLineRoot = nullptr;
    }

    /**
     * @brief Get the number of input sites
     */
    size_t VoronoiDiagram::getSiteCount() const
    {
        return sites.size();
    }

    /**
     * @brief Get the number of computed edges
     */
    size_t VoronoiDiagram::getEdgeCount() const
    {
        return edges.size();
    }

    // ============================================================================
    // PUBLIC METHODS - DIAGRAM ACCESS
    // ============================================================================

    /**
     * @brief Get the computed edges
     * @return Const reference to vector of edges
     */
    const std::vector<std::shared_ptr<Edge>>& VoronoiDiagram::getEdges() const
    {
        return edges;
    }

    /**
     * @brief Get the input sites
     * @return Const reference to vector of sites
     */
    const std::vector<Point>& VoronoiDiagram::getSites() const
    {
        return sites;
    }

    /**
     * @brief Set the bounding box for the diagram
     *
     * Infinite edges will be clipped to this bounding box.
     * Should be set before calling compute().
     *
     * @param minX Left boundary
     * @param maxX Right boundary
     * @param minY Bottom boundary
     * @param maxY Top boundary
     */
    void VoronoiDiagram::setBoundingBox(double minX, double maxX, double minY, double maxY)
    {
        // Store visible bounds (what the user sees)
        visibleMinX = minX;
        visibleMaxX = maxX;
        visibleMinY = minY;
        visibleMaxY = maxY;

        // Expand bounds for computation (allows edges to extend naturally)
        boundMinX = minX - COMPUTE_MARGIN;
        boundMaxX = maxX + COMPUTE_MARGIN;
        boundMinY = minY - COMPUTE_MARGIN;
        boundMaxY = maxY + COMPUTE_MARGIN;
    }

    // ============================================================================
    // MAIN ALGORITHM - COMPUTE
    // ============================================================================

    /**
     * @brief Main entry point for Fortune's algorithm
     *
     * This method implements the complete sweep line algorithm:
     * 1. Initialize by adding all sites as events
     * 2. Process events in order (highest y first)
     * 3. Handle site events by inserting new arcs
     * 4. Handle circle events by removing converging arcs
     * 5. Finish by clipping infinite edges to bounding box
     */
    void VoronoiDiagram::compute()
    {
        // Need at least 2 sites to have edges
        if (sites.size() < 2)
        {
            return;
        }

        // Clean up any previous computation
        cleanup();
        edges.clear();

        // Step 1: Initialize - add all sites to event queue
        initialize();

        // Step 2: Process all events
        while (!eventQueue.empty())
        {
            // Get the next event (highest y-coordinate)
            Event* event = eventQueue.top();
            eventQueue.pop();

            // Update sweep line position
            sweepLineY = event->point.y;

            // Skip invalidated circle events
            if (event->type == EventType::CIRCLE_EVENT && !event->isValid)
            {
                delete event;
                continue;
            }

            // Process the event based on its type
            if (event->type == EventType::SITE_EVENT)
            {
                handleSiteEvent(event);
            }
            else
            {
                handleCircleEvent(event);
            }

            delete event;
        }

        // Step 3: Finish incomplete edges by clipping to bounding box
        finishEdges();
    }

    // ============================================================================
    // INITIALIZATION
    // ============================================================================

    /**
     * @brief Initialize the algorithm by adding all site events
     *
     * Creates an event for each input site and adds it to the priority queue.
     * Events are automatically sorted by y-coordinate (highest first).
     */
    void VoronoiDiagram::initialize()
    {
        // Sort sites for consistent processing
        std::sort(sites.begin(), sites.end());

        // Add each site as a site event
        for (const Point& site : sites)
        {
            Event* event = new Event(site);
            eventQueue.push(event);
        }
    }

    // ============================================================================
    // EVENT HANDLING
    // ============================================================================

    /**
     * @brief Handle a site event
     *
     * When the sweep line encounters a new site:
     * 1. Find the arc directly above the new site
     * 2. Split that arc and insert a new arc for the new site
     * 3. Create edges between the new site and its neighbors
     * 4. Check for potential circle events
     *
     * @param event The site event to process
     */
    void VoronoiDiagram::handleSiteEvent(Event* event)
    {
        Point newSite = event->point;

        // Special case: first site - just create an arc
        if (beachLineRoot == nullptr)
        {
            beachLineRoot = new Arc(newSite);
            return;
        }

        // Find the arc directly above the new site
        Arc* arcAbove = findArcAbove(newSite.x);

        // Invalidate any circle event for the arc we're splitting
        invalidateCircleEvent(arcAbove);

        // Insert the new arc (splits arcAbove into three arcs)
        insertArc(newSite, arcAbove);
    }

    /**
     * @brief Handle a circle event
     *
     * When three consecutive arcs converge at a point:
     * 1. The middle arc disappears
     * 2. A Voronoi vertex is created at the convergence point
     * 3. The edges from the neighboring arcs meet at this vertex
     * 4. A new edge begins tracing between the remaining neighbors
     *
     * @param event The circle event to process
     */
    void VoronoiDiagram::handleCircleEvent(Event* event)
    {
        Arc* arc = event->arc;

        // Skip if this event was invalidated
        if (!event->isValid)
        {
            return;
        }

        // Get the vertex location (center of the circle)
        Point vertex;
        if (!getCircleCenter(arc->prev->site, arc->site, arc->next->site, vertex))
        {
            return;  // Shouldn't happen if event was valid
        }

        // Save neighbors BEFORE removing the arc
        Arc* prevArc = arc->prev;
        Arc* nextArc = arc->next;

        // Complete the edges that meet at this vertex
        if (arc->leftEdge)
        {
            arc->leftEdge->setEnd(vertex);
        }
        if (arc->rightEdge)
        {
            arc->rightEdge->setEnd(vertex);
        }

        // Create a new edge between the two remaining neighbors
        std::shared_ptr<Edge> newEdge = createEdge(prevArc->site, nextArc->site);
        newEdge->setStart(vertex);

        // Update edge pointers for remaining arcs
        prevArc->rightEdge = newEdge;
        nextArc->leftEdge = newEdge;

        // Invalidate circle events for neighboring arcs
        invalidateCircleEvent(prevArc);
        invalidateCircleEvent(nextArc);

        // Remove the disappearing arc (arc pointer becomes invalid after this)
        removeArc(arc);

        // Check for new circle events with the newly adjacent arcs
        checkCircleEvent(prevArc);
        checkCircleEvent(nextArc);
    }

    // ============================================================================
    // BEACH LINE MANIPULATION
    // ============================================================================

    /**
     * @brief Find the arc directly above a given x-coordinate
     *
     * Traverses the beach line to find which parabolic arc is directly
     * above the given x-coordinate at the current sweep line position.
     *
     * @param x The x-coordinate to search for
     * @return Pointer to the arc above x
     */
    Arc* VoronoiDiagram::findArcAbove(double x)
    {
        Arc* current = beachLineRoot;

        while (current != nullptr)
        {
            // Check left boundary
            if (current->prev != nullptr)
            {
                double leftX = getParabolaIntersectionX(current->prev->site, current->site, sweepLineY);
                if (x < leftX)
                {
                    current = current->prev;
                    continue;
                }
            }

            // Check right boundary
            if (current->next != nullptr)
            {
                double rightX = getParabolaIntersectionX(current->site, current->next->site, sweepLineY);
                if (x > rightX)
                {
                    current = current->next;
                    continue;
                }
            }

            // Found the arc
            break;
        }

        return current;
    }

    /**
     * @brief Insert a new arc into the beach line
     *
     * When a new site is encountered, we split the arc above it into three:
     * - Left portion of original arc
     * - New arc for the new site
     * - Right portion of original arc (copy of left portion's site)
     *
     * @param newSite The site for the new arc
     * @param arcAbove The arc being split
     */
    void VoronoiDiagram::insertArc(const Point& newSite, Arc* arcAbove)
    {
        // Create the three arcs
        Arc* leftArc = new Arc(arcAbove->site);    // Left portion
        Arc* middleArc = new Arc(newSite);         // New site's arc
        Arc* rightArc = new Arc(arcAbove->site);   // Right portion

        // Link them together
        leftArc->prev = arcAbove->prev;
        leftArc->next = middleArc;
        middleArc->prev = leftArc;
        middleArc->next = rightArc;
        rightArc->prev = middleArc;
        rightArc->next = arcAbove->next;

        // Update neighbors' pointers
        if (arcAbove->prev != nullptr)
        {
            arcAbove->prev->next = leftArc;
        }
        else
        {
            beachLineRoot = leftArc;  // New root
        }

        if (arcAbove->next != nullptr)
        {
            arcAbove->next->prev = rightArc;
        }

        // Calculate the start point for the new edges
        // This is where the new site's parabola first touches the arc above
        Point edgeStart;
        edgeStart.x = newSite.x;

        // Calculate y on the parabola of arcAbove at x = newSite.x
        double dy = arcAbove->site.y - sweepLineY;
        if (std::abs(dy) > EPSILON)
        {
            double dx = newSite.x - arcAbove->site.x;
            edgeStart.y = (dx * dx) / (2.0 * dy) + (arcAbove->site.y + sweepLineY) / 2.0;
        }
        else
        {
            // Degenerate case: site is on sweep line
            edgeStart.y = newSite.y;
        }

        // Create edges between the new site and the split site
        std::shared_ptr<Edge> leftEdge = createEdge(arcAbove->site, newSite);
        std::shared_ptr<Edge> rightEdge = createEdge(newSite, arcAbove->site);

        // Set the start points for both edges (they start at the same point)
        leftEdge->setStart(edgeStart);
        rightEdge->setStart(edgeStart);

        // Assign edges to arcs
        leftArc->rightEdge = leftEdge;
        middleArc->leftEdge = leftEdge;
        middleArc->rightEdge = rightEdge;
        rightArc->leftEdge = rightEdge;

        // Preserve edge pointers from original arc
        leftArc->leftEdge = arcAbove->leftEdge;
        rightArc->rightEdge = arcAbove->rightEdge;

        // Delete the original arc
        delete arcAbove;

        // Check for circle events
        checkCircleEvent(leftArc);
        checkCircleEvent(rightArc);

    }

    /**
     * @brief Remove an arc from the beach line
     *
     * Updates the doubly-linked list to remove the given arc.
     *
     * @param arc The arc to remove
     */
    void VoronoiDiagram::removeArc(Arc* arc)
    {
        // Update neighbors' pointers
        if (arc->prev != nullptr)
        {
            arc->prev->next = arc->next;
        }
        else
        {
            beachLineRoot = arc->next;
        }

        if (arc->next != nullptr)
        {
            arc->next->prev = arc->prev;
        }

        delete arc;
    }

    // ============================================================================
    // CIRCLE EVENT MANAGEMENT
    // ============================================================================

    /**
     * @brief Check if three consecutive arcs will create a circle event
     *
     * A circle event occurs when three consecutive arcs (left, middle, right)
     * are defined by sites that lie on a circle whose bottom touches the
     * sweep line at some point below the current position.
     *
     * @param arc The middle arc to check (with its neighbors)
     */
    void VoronoiDiagram::checkCircleEvent(Arc* arc)
    {
        // Need three consecutive arcs
        if (arc == nullptr || arc->prev == nullptr || arc->next == nullptr)
        {
            return;
        }

        Point a = arc->prev->site;
        Point b = arc->site;
        Point c = arc->next->site;

        // Check if the sites would create a valid circle event
        // The arcs must converge (right-hand turn test)
        double crossProduct = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        if (crossProduct >= 0)
        {
            return;  // Sites are collinear or turn left - no convergence
        }

        // Calculate the circle event position (bottom of circumcircle)
        double eventY = getCircleEventY(a, b, c);

        // Only create event if it's below current sweep line
        if (eventY < sweepLineY - EPSILON)
        {
            Point center;
            if (getCircleCenter(a, b, c, center))
            {
                // Event point is at the bottom of the circle
                Point eventPoint(center.x, eventY);
                Event* event = new Event(eventPoint, arc);
                arc->circleEvent = event;
                eventQueue.push(event);
            }
        }
    }

    /**
     * @brief Invalidate a circle event
     *
     * When an arc is modified, any pending circle event for that arc
     * becomes invalid and should not be processed.
     *
     * @param arc The arc whose circle event should be invalidated
     */
    void VoronoiDiagram::invalidateCircleEvent(Arc* arc)
    {
        if (arc != nullptr && arc->circleEvent != nullptr)
        {
            arc->circleEvent->isValid = false;
            arc->circleEvent = nullptr;
        }
    }

    // ============================================================================
    // GEOMETRY UTILITIES
    // ============================================================================

    /**
     * @brief Calculate x-coordinate where two parabolas intersect
     *
     * Each parabola is defined by a focus (site) and directrix (sweep line).
     * This function finds where two such parabolas meet.
     *
     * Mathematical derivation:
     * For a parabola with focus (fx, fy) and directrix y = ly:
     * Any point (x, y) on the parabola satisfies: (y - fy)^2 + (x - fx)^2 = (y - ly)^2
     * Solving for y: y = ((x - fx)^2 + fy^2 - ly^2) / (2 * (fy - ly))
     *
     * Setting two parabolas equal and solving for x gives the intersection.
     *
     * @param p1 Focus of first parabola
     * @param p2 Focus of second parabola
     * @param ly Y-coordinate of directrix
     * @return X-coordinate of intersection
     */
    double VoronoiDiagram::getParabolaIntersectionX(const Point& p1, const Point& p2, double ly)
    {
        // Handle degenerate cases
        if (isEqual(p1.y, p2.y))
        {
            // Both foci at same y - intersection is at midpoint
            return (p1.x + p2.x) / 2.0;
        }

        if (isEqual(p1.y, ly))
        {
            // p1 is on the sweep line - intersection is directly above p1
            return p1.x;
        }

        if (isEqual(p2.y, ly))
        {
            // p2 is on the sweep line - intersection is directly above p2
            return p2.x;
        }

        // General case - solve quadratic equation
        double d1 = 1.0 / (2.0 * (p1.y - ly));
        double d2 = 1.0 / (2.0 * (p2.y - ly));

        double a = d1 - d2;
        double b = 2.0 * (p2.x * d2 - p1.x * d1);
        double c = (p1.x * p1.x + p1.y * p1.y - ly * ly) * d1
            - (p2.x * p2.x + p2.y * p2.y - ly * ly) * d2;

        double discriminant = b * b - 4.0 * a * c;

        if (discriminant < 0)
        {
            discriminant = 0;  // Handle numerical errors
        }

        double x1 = (-b + std::sqrt(discriminant)) / (2.0 * a);
        double x2 = (-b - std::sqrt(discriminant)) / (2.0 * a);

        // Return the intersection point that's relevant
        // (the one between the two parabolas at this sweep line position)
        if (p1.y < p2.y)
        {
            return std::max(x1, x2);
        }
        else
        {
            return std::min(x1, x2);
        }
    }

    /**
     * @brief Calculate the center of a circle through three points
     *
     * Uses the circumcircle formula to find the center of the circle
     * passing through all three input points.
     *
     * @param a First point
     * @param b Second point
     * @param c Third point
     * @param center Output: center of the circle
     * @return True if circle exists (points not collinear)
     */
    bool VoronoiDiagram::getCircleCenter(const Point& a, const Point& b, const Point& c, Point& center)
    {
        double ax = a.x, ay = a.y;
        double bx = b.x, by = b.y;
        double cx = c.x, cy = c.y;

        double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

        if (isZero(d))
        {
            return false;  // Points are collinear
        }

        double aSq = ax * ax + ay * ay;
        double bSq = bx * bx + by * by;
        double cSq = cx * cx + cy * cy;

        center.x = (aSq * (by - cy) + bSq * (cy - ay) + cSq * (ay - by)) / d;
        center.y = (aSq * (cx - bx) + bSq * (ax - cx) + cSq * (bx - ax)) / d;

        return true;
    }

    /**
     * @brief Get the y-coordinate where a circle event will occur
     *
     * The circle event occurs when the sweep line is tangent to the
     * bottom of the circumcircle of the three sites.
     *
     * @param a First site
     * @param b Second site (middle)
     * @param c Third site
     * @return Y-coordinate of the circle event
     */
    double VoronoiDiagram::getCircleEventY(const Point& a, const Point& b, const Point& c)
    {
        Point center;
        if (!getCircleCenter(a, b, c, center))
        {
            return -INF;  // No valid circle
        }

        // Radius of the circle
        double radius = center.distanceTo(a);

        // Circle event occurs at bottom of circle
        return center.y - radius;
    }

    // ============================================================================
    // EDGE MANAGEMENT
    // ============================================================================

    /**
     * @brief Create a new edge between two sites
     *
     * @param site1 First site
     * @param site2 Second site
     * @return Shared pointer to the new edge
     */
    std::shared_ptr<Edge> VoronoiDiagram::createEdge(const Point& site1, const Point& site2)
    {
        auto edge = std::make_shared<Edge>(site1, site2);
        edges.push_back(edge);
        return edge;
    }

    /**
     * @brief Finish all incomplete edges
     *
     * After processing all events, some edges may still be incomplete
     * (extending to infinity). This method clips them to the bounding box.
     */
    void VoronoiDiagram::finishEdges()
    {
        // Step 1: Extend all edges that need it (incomplete OR have invalid end point)
        for (auto& edge : edges)
        {
            bool hasValidEnd = (std::abs(edge->end.x) > EPSILON || std::abs(edge->end.y) > EPSILON);
            if (!edge->isComplete || !hasValidEnd)
            {
                extendEdgeToComputeBounds(edge);
            }
        }

        // Step 2: Clip all edges to the VISIBLE bounds
        for (auto& edge : edges)
        {
            clipEdgeToVisibleBounds(edge);
        }

        // Step 3: Remove invalid edges
        edges.erase(
            std::remove_if(edges.begin(), edges.end(),
                           [](const std::shared_ptr<Edge>& e) {
                               double len = std::sqrt(std::pow(e->end.x - e->start.x, 2) +
                                                      std::pow(e->end.y - e->start.y, 2));
                               if (len < EPSILON) return true;

                               if (std::abs(e->start.x) < EPSILON && std::abs(e->start.y) < EPSILON &&
                                   std::abs(e->end.x) < EPSILON && std::abs(e->end.y) < EPSILON)
                               {
                                   return true;
                               }

                               return false;
                           }),
            edges.end()
        );
    }

    /**
     * @brief Extend an incomplete edge to the computation bounding box
     *
     * For edges that have a start point but no end point, this extends
     * them in the appropriate direction until they hit the compute bounds.
     *
     * @param edge The edge to extend
     */
    void VoronoiDiagram::extendEdgeToComputeBounds(std::shared_ptr<Edge>& edge)
    {
        // Check if edge has a valid start point
        bool hasValidStart = (std::abs(edge->start.x) > EPSILON || std::abs(edge->start.y) > EPSILON);

        if (!hasValidStart)
        {
            // No valid start point - mark as invalid
            edge->start = Point(0, 0);
            edge->end = Point(0, 0);
            edge->isComplete = true;
            return;
        }

        // Direction perpendicular to the line between sites
        double dx = edge->site2.y - edge->site1.y;
        double dy = -(edge->site2.x - edge->site1.x);

        // Normalize direction
        double len = std::sqrt(dx * dx + dy * dy);
        if (len < EPSILON)
        {
            edge->start = Point(0, 0);
            edge->end = Point(0, 0);
            edge->isComplete = true;
            return;
        }
        dx /= len;
        dy /= len;

        // Find intersection with COMPUTE bounding box
        double tMax = INF;

        if (std::abs(dx) > EPSILON)
        {
            double t1 = (boundMinX - edge->start.x) / dx;
            double t2 = (boundMaxX - edge->start.x) / dx;
            if (t1 > EPSILON) tMax = std::min(tMax, t1);
            if (t2 > EPSILON) tMax = std::min(tMax, t2);
        }

        if (std::abs(dy) > EPSILON)
        {
            double t1 = (boundMinY - edge->start.y) / dy;
            double t2 = (boundMaxY - edge->start.y) / dy;
            if (t1 > EPSILON) tMax = std::min(tMax, t1);
            if (t2 > EPSILON) tMax = std::min(tMax, t2);
        }

        if (tMax == INF || tMax < EPSILON)
        {
            edge->start = Point(0, 0);
            edge->end = Point(0, 0);
            edge->isComplete = true;
            return;
        }

        // Calculate endpoint on compute bounds
        double endX = edge->start.x + dx * tMax;
        double endY = edge->start.y + dy * tMax;

        edge->end = Point(endX, endY);
        edge->isComplete = true;
    }

    /**
     * @brief Clip an edge to the visible bounding box
     *
     * Uses Cohen-Sutherland algorithm to clip the edge to what the user sees.
     *
     * @param edge The edge to clip
     */
    void VoronoiDiagram::clipEdgeToVisibleBounds(std::shared_ptr<Edge>& edge)
    {
        // Check for invalid endpoints (0,0 is our "unset" marker)
        bool startValid = (std::abs(edge->start.x) > EPSILON || std::abs(edge->start.y) > EPSILON);
        bool endValid = (std::abs(edge->end.x) > EPSILON || std::abs(edge->end.y) > EPSILON);

        if (!startValid || !endValid)
        {
            // Edge has invalid endpoint - mark for removal
            edge->start = Point(0, 0);
            edge->end = Point(0, 0);
            return;
        }

        double x0 = edge->start.x;
        double y0 = edge->start.y;
        double x1 = edge->end.x;
        double y1 = edge->end.y;

        // Cohen-Sutherland outcodes (using VISIBLE bounds)
        const int INSIDE = 0;
        const int LEFT = 1;
        const int RIGHT = 2;
        const int BOTTOM = 4;
        const int TOP = 8;

        auto computeOutCode = [this](double x, double y) {
            int code = INSIDE;
            if (x < visibleMinX) code |= LEFT;
            else if (x > visibleMaxX) code |= RIGHT;
            if (y < visibleMinY) code |= TOP;
            else if (y > visibleMaxY) code |= BOTTOM;
            return code;
            };

        int outcode0 = computeOutCode(x0, y0);
        int outcode1 = computeOutCode(x1, y1);

        while (true)
        {
            if (!(outcode0 | outcode1))
            {
                // Both inside - accept
                break;
            }
            else if (outcode0 & outcode1)
            {
                // Both outside same region - reject (mark as invalid)
                edge->start = Point(0, 0);
                edge->end = Point(0, 0);
                return;
            }
            else
            {
                // Needs clipping
                int outcodeOut = outcode0 ? outcode0 : outcode1;
                double x, y;

                if (outcodeOut & BOTTOM)
                {
                    x = x0 + (x1 - x0) * (visibleMaxY - y0) / (y1 - y0);
                    y = visibleMaxY;
                }
                else if (outcodeOut & TOP)
                {
                    x = x0 + (x1 - x0) * (visibleMinY - y0) / (y1 - y0);
                    y = visibleMinY;
                }
                else if (outcodeOut & RIGHT)
                {
                    y = y0 + (y1 - y0) * (visibleMaxX - x0) / (x1 - x0);
                    x = visibleMaxX;
                }
                else
                { // LEFT
                    y = y0 + (y1 - y0) * (visibleMinX - x0) / (x1 - x0);
                    x = visibleMinX;
                }

                if (outcodeOut == outcode0)
                {
                    x0 = x;
                    y0 = y;
                    outcode0 = computeOutCode(x0, y0);
                }
                else
                {
                    x1 = x;
                    y1 = y;
                    outcode1 = computeOutCode(x1, y1);
                }
            }
        }

        edge->start = Point(x0, y0);
        edge->end = Point(x1, y1);
    }

    /**
     * @brief Clip a line segment to a rectangular bounding box using Cohen-Sutherland algorithm
     *
     * @param x0, y0 Start point (modified in place)
     * @param x1, y1 End point (modified in place)
     * @return true if any part of line is inside box, false if entirely outside
     */
    bool VoronoiDiagram::clipLineToBox(double& x0, double& y0, double& x1, double& y1)
    {
        // Cohen-Sutherland outcodes
        const int INSIDE = 0;
        const int LEFT = 1;
        const int RIGHT = 2;
        const int BOTTOM = 4;
        const int TOP = 8;

        auto computeOutCode = [this](double x, double y) {
            int code = INSIDE;
            if (x < boundMinX) code |= LEFT;
            else if (x > boundMaxX) code |= RIGHT;
            if (y < boundMinY) code |= TOP;      // Note: y increases downward in screen coords
            else if (y > boundMaxY) code |= BOTTOM;
            return code;
            };

        int outcode0 = computeOutCode(x0, y0);
        int outcode1 = computeOutCode(x1, y1);

        while (true)
        {
            if (!(outcode0 | outcode1))
            {
                // Both inside
                return true;
            }
            else if (outcode0 & outcode1)
            {
                // Both outside same region
                return false;
            }
            else
            {
                // Needs clipping
                int outcodeOut = outcode0 ? outcode0 : outcode1;
                double x, y;

                if (outcodeOut & BOTTOM)
                {
                    x = x0 + (x1 - x0) * (boundMaxY - y0) / (y1 - y0);
                    y = boundMaxY;
                }
                else if (outcodeOut & TOP)
                {
                    x = x0 + (x1 - x0) * (boundMinY - y0) / (y1 - y0);
                    y = boundMinY;
                }
                else if (outcodeOut & RIGHT)
                {
                    y = y0 + (y1 - y0) * (boundMaxX - x0) / (x1 - x0);
                    x = boundMaxX;
                }
                else
                { // LEFT
                    y = y0 + (y1 - y0) * (boundMinX - x0) / (x1 - x0);
                    x = boundMinX;
                }

                if (outcodeOut == outcode0)
                {
                    x0 = x;
                    y0 = y;
                    outcode0 = computeOutCode(x0, y0);
                }
                else
                {
                    x1 = x;
                    y1 = y;
                    outcode1 = computeOutCode(x1, y1);
                }
            }
        }
    }

    // ============================================================================
    // FILE I/O
    // ============================================================================

    /**
     * @brief Export the diagram to a JSON file
     *
     * Creates a JSON file containing all sites and edges for visualization.
     *
     * @param filename Output filename
     * @return True if successful
     */
    bool VoronoiDiagram::exportToFile(const std::string& filename) const
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
            return false;
        }

        file << "{\n";

        // Export bounding box
        file << "  \"bounds\": {\n";
        file << "    \"minX\": " << visibleMinX << ",\n";
        file << "    \"maxX\": " << visibleMaxX << ",\n";
        file << "    \"minY\": " << visibleMinY << ",\n";
        file << "    \"maxY\": " << visibleMaxY << "\n";
        file << "  },\n";

        // Export sites
        file << "  \"sites\": [\n";
        for (size_t i = 0; i < sites.size(); ++i)
        {
            file << "    {\"x\": " << sites[i].x
                << ", \"y\": " << sites[i].y
                << ", \"id\": " << sites[i].id << "}";
            if (i < sites.size() - 1) file << ",";
            file << "\n";
        }
        file << "  ],\n";

        // Export edges
        file << "  \"edges\": [\n";
        for (size_t i = 0; i < edges.size(); ++i)
        {
            const auto& edge = edges[i];
            file << "    {\"start\": {\"x\": " << edge->start.x << ", \"y\": " << edge->start.y << "}, "
                << "\"end\": {\"x\": " << edge->end.x << ", \"y\": " << edge->end.y << "}, "
                << "\"site1_id\": " << edge->site1.id << ", "
                << "\"site2_id\": " << edge->site2.id << "}";
            if (i < edges.size() - 1) file << ",";
            file << "\n";
        }
        file << "  ]\n";

        file << "}\n";
        file.close();

        return true;
    }

    /**
     * @brief Import sites from a text file
     *
     * Reads sites from a file with one point per line (format: x y)
     *
     * @param filename Input filename
     * @return True if successful
     */
    bool VoronoiDiagram::importFromFile(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Could not open file for reading: " << filename << std::endl;
            return false;
        }

        clear();

        std::string line;
        while (std::getline(file, line))
        {
            // Skip empty lines and comments
            if (line.empty() || line[0] == '#')
            {
                continue;
            }

            std::istringstream iss(line);
            double x, y;
            if (iss >> x >> y)
            {
                addSite(x, y);
            }
        }

        file.close();
        return true;
    }

    // ============================================================================
    // CLEANUP
    // ============================================================================

    /**
     * @brief Free all dynamically allocated memory
     *
     * Cleans up the beach line and any remaining events in the queue.
     */
    void VoronoiDiagram::cleanup()
    {
        // Clean up beach line
        Arc* current = beachLineRoot;
        while (current != nullptr)
        {
            Arc* next = current->next;
            delete current;
            current = next;
        }
        beachLineRoot = nullptr;

        // Clean up event queue
        while (!eventQueue.empty())
        {
            Event* event = eventQueue.top();
            eventQueue.pop();
            delete event;
        }
    }

} // namespace Voronoi
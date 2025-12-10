/**
 * @file voronoi.h
 * @brief Header file for Fortune's Algorithm Voronoi Diagram Implementation
 *
 * This file contains all the data structures and class declarations needed
 * to implement Fortune's sweep line algorithm for computing Voronoi diagrams.
 *
 * Fortune's Algorithm Overview:
 * - Sweeps a horizontal line from top to bottom across the plane
 * - Maintains a "beach line" of parabolas representing equidistant points
 * - Processes two types of events: Site events and Circle events
 * - Produces a Voronoi diagram in O(n log n) time
 *
 * @author Jack Woods, Ashleigh Kirkpatrick
 * @date 2025
 */

#ifndef VORONOI_H
#define VORONOI_H

#include <vector>
#include <queue>
#include <set>
#include <memory>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace Voronoi
{

    // ============================================================================
    // FORWARD DECLARATIONS
    // ============================================================================

    struct Point;
    struct Edge;
    struct Arc;
    struct Event;
    class VoronoiDiagram;

    // ============================================================================
    // CONSTANTS
    // ============================================================================

    // Small epsilon value for floating point comparisons
    constexpr double EPSILON = 1e-9;

    // Large value representing infinity for unbounded edges
    constexpr double INF = 1e10;

    // ============================================================================
    // BASIC DATA STRUCTURES
    // ============================================================================

    /**
     * @struct Point
     * @brief Represents a 2D point with x and y coordinates
     *
     * Used for both input sites and vertices of the Voronoi diagram.
     * Includes basic comparison operators for use in sorted containers.
     */
    struct Point
    {
        double x;  // X coordinate
        double y;  // Y coordinate
        int id;    // Unique identifier for the site (used for output)

        /**
         * @brief Default constructor - initializes to origin
         */
        Point() : x(0), y(0), id(-1)
        {
        }

        /**
         * @brief Parameterized constructor
         * @param x X coordinate
         * @param y Y coordinate
         * @param id Optional identifier for the point
         */
        Point(double x, double y, int id = -1) : x(x), y(y), id(id)
        {
        }

        /**
         * @brief Equality comparison (within epsilon tolerance)
         */
        bool operator==(const Point& other) const
        {
            return std::abs(x - other.x) < EPSILON &&
                std::abs(y - other.y) < EPSILON;
        }

        /**
         * @brief Inequality comparison
         */
        bool operator!=(const Point& other) const
        {
            return !(*this == other);
        }

        /**
         * @brief Less than comparison (for sorting: by y, then by x)
         * Note: We sort by y descending because sweep line moves top to bottom
         */
        bool operator<(const Point& other) const
        {
            if (std::abs(y - other.y) > EPSILON)
            {
                return y > other.y;  // Higher y values come first
            }
            return x < other.x;
        }

        /**
         * @brief Calculate Euclidean distance to another point
         */
        double distanceTo(const Point& other) const
        {
            double dx = x - other.x;
            double dy = y - other.y;
            return std::sqrt(dx * dx + dy * dy);
        }
    };

    /**
     * @struct Edge
     * @brief Represents an edge in the Voronoi diagram
     *
     * Each edge is a bisector between two adjacent sites. Edges may be
     * half-infinite (extending to the bounding box) or finite.
     */
    struct Edge
    {
        Point start;       // Starting vertex of the edge
        Point end;         // Ending vertex of the edge
        Point site1;       // One of the two sites this edge bisects
        Point site2;       // The other site this edge bisects
        bool isComplete;   // True if both endpoints are known

        /**
         * @brief Default constructor
         */
        Edge() : isComplete(false)
        {
        }

        /**
         * @brief Constructor with sites
         * @param s1 First site
         * @param s2 Second site
         */
        Edge(const Point& s1, const Point& s2)
            : site1(s1), site2(s2), isComplete(false)
        {
        }

        /**
         * @brief Set the start point of the edge
         */
        void setStart(const Point& p)
        {
            start = p;
        }

        /**
         * @brief Set the end point and mark edge as complete
         */
        void setEnd(const Point& p)
        {
            end = p;
            isComplete = true;
        }
    };

    // ============================================================================
    // BEACH LINE STRUCTURES
    // ============================================================================

    /**
     * @struct Arc
     * @brief Represents a parabolic arc in the beach line
     *
     * The beach line is a sequence of parabolic arcs. Each arc is defined by
     * a site (focus) and the current position of the sweep line (directrix).
     * Arcs are stored in a doubly-linked list for efficient traversal.
     */
    struct Arc
    {
        Point site;                      // The site (focus) defining this parabola
        Arc* prev;                       // Previous arc in beach line (left neighbor)
        Arc* next;                       // Next arc in beach line (right neighbor)
        Event* circleEvent;              // Pointer to this arc's circle event (if any)
        std::shared_ptr<Edge> leftEdge;  // Edge being traced on the left
        std::shared_ptr<Edge> rightEdge; // Edge being traced on the right

        /**
         * @brief Constructor
         * @param s The site that defines this parabolic arc
         */
        Arc(const Point& s)
            : site(s), prev(nullptr), next(nullptr), circleEvent(nullptr)
        {
        }
    };

    // ============================================================================
    // EVENT STRUCTURES
    // ============================================================================

    /**
     * @enum EventType
     * @brief Types of events in Fortune's algorithm
     */
    enum class EventType
    {
        SITE_EVENT,    // A new site is encountered by the sweep line
        CIRCLE_EVENT   // Three consecutive arcs converge at a point
    };

    /**
     * @struct Event
     * @brief Represents an event in the priority queue
     *
     * Site events occur when the sweep line hits a new site.
     * Circle events occur when three consecutive arcs converge.
     */
    struct Event
    {
        EventType type;      // Type of this event
        Point point;         // Location where event occurs
        Arc* arc;            // Associated arc (for circle events)
        bool isValid;        // False if this event has been invalidated

        /**
         * @brief Constructor for site events
         */
        Event(const Point& p)
            : type(EventType::SITE_EVENT), point(p), arc(nullptr), isValid(true)
        {
        }

        /**
         * @brief Constructor for circle events
         */
        Event(const Point& p, Arc* a)
            : type(EventType::CIRCLE_EVENT), point(p), arc(a), isValid(true)
        {
        }

        /**
         * @brief Comparison for priority queue (process highest y first)
         */
        bool operator>(const Event& other) const
        {
            if (std::abs(point.y - other.point.y) > EPSILON)
            {
                return point.y < other.point.y;  // Lower y = lower priority
            }
            return point.x > other.point.x;
        }
    };

    /**
     * @struct EventComparator
     * @brief Comparator for the event priority queue
     */
    struct EventComparator
    {
        bool operator()(Event* a, Event* b) const
        {
            return *a > *b;
        }
    };

    // ============================================================================
    // MAIN VORONOI CLASS
    // ============================================================================

    /**
     * @class VoronoiDiagram
     * @brief Main class implementing Fortune's algorithm
     *
     * This class takes a set of input points (sites) and computes the Voronoi
     * diagram using Fortune's sweep line algorithm. The result is a set of
     * edges that partition the plane into regions closest to each site.
     */
    class VoronoiDiagram
    {
    public:
        // ========================================================================
        // PUBLIC METHODS
        // ========================================================================

        /**
         * @brief Default constructor
         */
        VoronoiDiagram();

        /**
         * @brief Destructor - cleans up allocated memory
         */
        ~VoronoiDiagram();

        /**
         * @brief Add a site to the diagram
         * @param x X coordinate of the site
         * @param y Y coordinate of the site
         */
        void addSite(double x, double y);

        /**
         * @brief Add a site using a Point object
         * @param p The point to add as a site
         */
        void addSite(const Point& p);

        /**
         * @brief Clear all sites and computed data
         */
        void clear();

        /**
         * @brief Compute the Voronoi diagram
         *
         * Main algorithm entry point. Processes all events in the priority
         * queue and builds the diagram edge by edge.
         */
        void compute();

        /**
         * @brief Get the computed edges
         * @return Vector of shared pointers to edges
         */
        const std::vector<std::shared_ptr<Edge>>& getEdges() const;

        /**
         * @brief Get the input sites
         * @return Vector of site points
         */
        const std::vector<Point>& getSites() const;

        /**
         * @brief Set the bounding box for clipping edges
         * @param minX Minimum x coordinate
         * @param maxX Maximum x coordinate
         * @param minY Minimum y coordinate
         * @param maxY Maximum y coordinate
         */
        void setBoundingBox(double minX, double maxX, double minY, double maxY);

        /**
         * @brief Export the diagram to a file (JSON format)
         * @param filename Output filename
         * @return True if export was successful
         */
        bool exportToFile(const std::string& filename) const;

        /**
         * @brief Import sites from a file
         * @param filename Input filename (one point per line: x y)
         * @return True if import was successful
         */
        bool importFromFile(const std::string& filename);

        /**
         * @brief Get the number of sites
         * @return Number of input sites
         */
        size_t getSiteCount() const;

        /**
         * @brief Get the number of edges
         * @return Number of computed edges
         */
        size_t getEdgeCount() const;

    private:
        // ========================================================================
        // PRIVATE DATA MEMBERS
        // ========================================================================

        std::vector<Point> sites;                    // Input sites
        std::vector<std::shared_ptr<Edge>> edges;    // Output edges
        Arc* beachLineRoot;                          // Root of beach line (doubly-linked list)

        // Priority queue for events (max-heap by y-coordinate)
        std::priority_queue<Event*, std::vector<Event*>, EventComparator> eventQueue;

        // Visible bounding box (what the user sees)
        double visibleMinX, visibleMaxX, visibleMinY, visibleMaxY;

        // Computation bounding box (expanded with margin for cleaner edge handling)
        double boundMinX, boundMaxX, boundMinY, boundMaxY;

        // Margin size for computation bounds
        static constexpr double COMPUTE_MARGIN = 200.0;

        // Current sweep line position
        double sweepLineY;

        // ========================================================================
        // PRIVATE METHODS - EVENT HANDLING
        // ========================================================================

        /**
         * @brief Initialize the algorithm by adding all site events
         */
        void initialize();

        /**
         * @brief Process a site event
         * @param event The site event to process
         */
        void handleSiteEvent(Event* event);

        /**
         * @brief Process a circle event
         * @param event The circle event to process
         */
        void handleCircleEvent(Event* event);

        // ========================================================================
        // PRIVATE METHODS - BEACH LINE MANIPULATION
        // ========================================================================

        /**
         * @brief Find the arc directly above a point
         * @param x X coordinate to search for
         * @return Pointer to the arc above x
         */
        Arc* findArcAbove(double x);

        /**
         * @brief Insert a new arc into the beach line
         * @param newSite The site for the new arc
         * @param arcAbove The arc being split
         */
        void insertArc(const Point& newSite, Arc* arcAbove);

        /**
         * @brief Remove an arc from the beach line
         * @param arc The arc to remove
         */
        void removeArc(Arc* arc);

        // ========================================================================
        // PRIVATE METHODS - CIRCLE EVENTS
        // ========================================================================

        /**
         * @brief Check if three consecutive arcs will converge
         * @param arc The middle arc to check
         */
        void checkCircleEvent(Arc* arc);

        /**
         * @brief Invalidate a circle event (mark as not to be processed)
         * @param arc The arc whose event should be invalidated
         */
        void invalidateCircleEvent(Arc* arc);

        // ========================================================================
        // PRIVATE METHODS - GEOMETRY UTILITIES
        // ========================================================================

        /**
         * @brief Compute the x-coordinate of intersection between two parabolas
         * @param p1 Focus of first parabola
         * @param p2 Focus of second parabola
         * @param ly Y-coordinate of the directrix (sweep line)
         * @return X-coordinate of intersection
         */
        double getParabolaIntersectionX(const Point& p1, const Point& p2, double ly);

        /**
         * @brief Compute the center of a circle through three points
         * @param a First point
         * @param b Second point
         * @param c Third point
         * @param center Output: center of the circle
         * @return True if circle exists (points not collinear)
         */
        bool getCircleCenter(const Point& a, const Point& b, const Point& c, Point& center);

        /**
         * @brief Get the y-coordinate of circle event for three points
         * @param a First point
         * @param b Second point (middle)
         * @param c Third point
         * @return Y-coordinate where circle touches sweep line (bottom of circle)
         */
        double getCircleEventY(const Point& a, const Point& b, const Point& c);

        // ========================================================================
        // PRIVATE METHODS - EDGE MANAGEMENT
        // ========================================================================

        /**
         * @brief Create a new edge between two sites
         * @param site1 First site
         * @param site2 Second site
         * @return Shared pointer to new edge
         */
        std::shared_ptr<Edge> createEdge(const Point& site1, const Point& site2);

        /**
         * @brief Finish all incomplete edges (clip to bounding box)
         */
        void finishEdges();

        /**
         * @brief Extend an incomplete edge to the computation bounding box
         * @param edge The edge to extend
         */
        void extendEdgeToComputeBounds(std::shared_ptr<Edge>& edge);

        /**
         * @brief Clip an edge to the visible bounding box
         * @param edge The edge to clip
         */
        void clipEdgeToVisibleBounds(std::shared_ptr<Edge>& edge);

        /**
         * @brief Clip a line segment to the bounding box (Cohen-Sutherland)
         * @param x0, y0, x1, y1 Line endpoints (modified in place)
         * @return True if line intersects bounding box
         */
        bool clipLineToBox(double& x0, double& y0, double& x1, double& y1);

        // ========================================================================
        // PRIVATE METHODS - CLEANUP
        // ========================================================================

        /**
         * @brief Free all allocated memory
         */
        void cleanup();
    };

    // ============================================================================
    // UTILITY FUNCTIONS
    // ============================================================================

    /**
     * @brief Check if a value is approximately zero
     */
    inline bool isZero(double val)
    {
        return std::abs(val) < EPSILON;
    }

    /**
     * @brief Check if two doubles are approximately equal
     */
    inline bool isEqual(double a, double b)
    {
        return std::abs(a - b) < EPSILON;
    }

    /**
     * @brief Clamp a value to a range
     */
    inline double clamp(double val, double minVal, double maxVal)
    {
        return std::max(minVal, std::min(val, maxVal));
    }

} // namespace Voronoi

#endif // VORONOI_H
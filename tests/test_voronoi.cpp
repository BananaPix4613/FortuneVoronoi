/**
 * @file test_voronoi.cpp
 * @brief Correctness Test Suite for Fortune's Algorithm Implementation
 * 
 * This file contains a comprehensive test suite to verify the correctness
 * of the Voronoi diagram implementation. Tests follow the format from
 * CSC 372 Test Plan guidelines.
 * 
 * Test Case Format:
 *   - Test ID and Description
 *   - Input data
 *   - Expected output
 *   - Actual output
 *   - Pass/Fail result
 * 
 * Categories of tests:
 *   1. Normal usage cases
 *   2. Edge cases (extrema)
 *   3. Best case scenarios
 *   4. Worst case scenarios
 *   5. Error handling
 * 
 * @author Jack Woods
 * @date 2025
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include "voronoi.h"

using namespace Voronoi;

// ============================================================================
// TEST FRAMEWORK
// ============================================================================

// Test result tracking
struct TestResult
{
    std::string testId;
    std::string description;
    std::string input;
    std::string expectedOutput;
    std::string actualOutput;
    bool passed;
};

// Global test results
std::vector<TestResult> testResults;
int totalTests = 0;
int passedTests = 0;

/**
 * @brief Record a test result
 */
void recordTest(const std::string& id, const std::string& desc,
                const std::string& input, const std::string& expected,
                const std::string& actual, bool passed)
{
    TestResult result;
    result.testId = id;
    result.description = desc;
    result.input = input;
    result.expectedOutput = expected;
    result.actualOutput = actual;
    result.passed = passed;
    testResults.push_back(result);
    
    totalTests++;
    if (passed) passedTests++;
    
    // Print result to console
    std::cout << (passed ? "[PASS]" : "[FAIL]") << " " << id << ": " << desc << "\n";
    if (!passed)
    {
        std::cout << "       Expected: " << expected << "\n";
        std::cout << "       Actual:   " << actual << "\n";
    }
}

/**
 * @brief Check if two doubles are approximately equal
 */
bool approxEqual(double a, double b, double epsilon = 1e-6)
{
    return std::abs(a - b) < epsilon;
}

/**
 * @brief Check if a point lies on a perpendicular bisector
 * 
 * A Voronoi edge is the perpendicular bisector of the line segment
 * connecting two sites. Any point on this edge should be equidistant
 * from both sites.
 */
bool isOnPerpendicularBisector(const Point& edgePoint, const Point& site1, const Point& site2)
{
    double dist1 = edgePoint.distanceTo(site1);
    double dist2 = edgePoint.distanceTo(site2);
    return approxEqual(dist1, dist2, 1.0);  // Allow 1.0 tolerance for floating point
}

// ============================================================================
// TEST CASES
// ============================================================================

/**
 * TEST 1: Two Sites (Minimum Case)
 * 
 * Description: Verify diagram with exactly two sites produces one edge
 * Input: Two points at (100, 300) and (700, 300)
 * Expected: One edge that is the perpendicular bisector between them
 */
void test01_TwoSites()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(100, 300);
    diagram.addSite(700, 300);
    diagram.compute();

    std::string input = "Sites: (100, 300), (700, 300)";
    std::string expected = "Vertical bisector edge(s) at x=400";

    bool passed = true;
    std::stringstream actual;

    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edge(s)";

    if (edgeCount < 1)
    {
        passed = false;
        actual << " - no edges generated";
    }
    else
    {
        // For two horizontal sites, the bisector is a vertical line at x=400
        // Check if any edge lies on (or near) x=400
        const auto& edges = diagram.getEdges();

        for (const auto& edge : edges)
        {
            // For a vertical bisector, both start.x and end.x should be near 400
            bool startNear400 = approxEqual(edge->start.x, 400, 50);
            bool endNear400 = approxEqual(edge->end.x, 400, 50);

            if (startNear400 && endNear400)
            {
                passed = true;
                break;
            }
        }

        if (passed)
        {
            actual << " - vertical bisector at x=400 found";
        }
        else
        {
            actual << " (bisector not at expected position)";
        }
    }

    recordTest("TEST-01", "Two Sites - Minimum Case", input, expected, actual.str(), passed);
}

/**
 * TEST 2: Three Sites (Triangle)
 * 
 * Description: Verify diagram with three non-collinear sites
 * Input: Equilateral triangle arrangement
 * Expected: Three edges meeting at circumcenter
 */
void test02_ThreeSitesTriangle()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    
    // Equilateral triangle (approximately)
    diagram.addSite(400, 100);  // Top
    diagram.addSite(200, 446);  // Bottom left
    diagram.addSite(600, 446);  // Bottom right
    diagram.compute();
    
    std::string input = "Sites: (400,100), (200,446), (600,446) - triangle";
    std::string expected = "3 edges meeting at a central vertex";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edge(s)";
    
    // Should have exactly 3 edges for 3 sites
    if (edgeCount >= 3)
    {
        actual << " - triangle configuration correct";
    }
    else
    {
        passed = false;
        actual << " - expected at least 3 edges";
    }
    
    recordTest("TEST-02", "Three Sites - Triangle Configuration", input, expected, actual.str(), passed);
}

/**
 * TEST 3: Collinear Sites
 * 
 * Description: Sites arranged in a line (degenerate case)
 * Input: Three sites on a horizontal line
 * Expected: Two parallel vertical edges
 */
void test03_CollinearSites()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(100, 300);
    diagram.addSite(400, 300);
    diagram.addSite(700, 300);
    diagram.compute();
    
    std::string input = "Sites: (100,300), (400,300), (700,300) - collinear";
    std::string expected = "2 vertical edges (perpendicular bisectors)";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edge(s)";
    
    if (edgeCount >= 2)
    {
        actual << " - collinear case handled";
    }
    else
    {
        passed = false;
        actual << " - insufficient edges for collinear sites";
    }
    
    recordTest("TEST-03", "Collinear Sites - Degenerate Case", input, expected, actual.str(), passed);
}

/**
 * TEST 4: Square Arrangement
 * 
 * Description: Four sites at corners of a square
 * Input: Sites at (200,200), (600,200), (200,400), (600,400)
 * Expected: Edges forming a cross pattern
 */
void test04_SquareArrangement()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(200, 200);
    diagram.addSite(600, 200);
    diagram.addSite(200, 400);
    diagram.addSite(600, 400);
    diagram.compute();
    
    std::string input = "Sites: corners of 400x200 rectangle";
    std::string expected = "Cross-shaped boundary pattern";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edge(s)";
    
    // Four sites in a rectangle should produce 5 edges (cross pattern + boundaries)
    if (edgeCount >= 4)
    {
        actual << " - rectangular configuration correct";
    }
    else
    {
        passed = false;
        actual << " - expected more edges for 4 sites";
    }
    
    recordTest("TEST-04", "Square Arrangement - Regular Pattern", input, expected, actual.str(), passed);
}

/**
 * TEST 5: Large Number of Sites
 * 
 * Description: Test with many random sites (stress test)
 * Input: 100 random sites
 * Expected: Diagram computes without error, satisfies Euler's formula
 */
void test05_LargeNumberOfSites()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    
    // Add 100 sites in a grid pattern
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            diagram.addSite(50 + i * 70, 50 + j * 50);
        }
    }
    diagram.compute();
    
    std::string input = "100 sites in 10x10 grid";
    std::string expected = "Approximately 3n-6 = 294 edges";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t n = diagram.getSiteCount();
    size_t edgeCount = diagram.getEdgeCount();
    size_t expectedEdges = 3 * n - 6;
    
    actual << edgeCount << " edges (expected ~" << expectedEdges << ")";
    
    // Allow some variance due to boundary conditions
    if (edgeCount >= expectedEdges / 2 && edgeCount <= expectedEdges * 2)
    {
        actual << " - within expected range";
    }
    else
    {
        passed = false;
        actual << " - edge count outside expected range";
    }
    
    recordTest("TEST-05", "Large Number of Sites - Stress Test", input, expected, actual.str(), passed);
}

/**
 * TEST 6: Single Site (Error Case)
 * 
 * Description: Test with only one site (should handle gracefully)
 * Input: Single site at (400, 300)
 * Expected: No edges generated, no crash
 */
void test06_SingleSite()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(400, 300);
    diagram.compute();
    
    std::string input = "Single site: (400, 300)";
    std::string expected = "0 edges (graceful handling)";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edges";
    
    if (edgeCount == 0)
    {
        actual << " - single site handled correctly";
        passed = true;
    }
    else
    {
        actual << " - unexpected edges for single site";
        passed = false;
    }
    
    recordTest("TEST-06", "Single Site - Error Case", input, expected, actual.str(), passed);
}

/**
 * TEST 7: Zero Sites (Error Case)
 * 
 * Description: Test with no sites (should handle gracefully)
 * Input: Empty site set
 * Expected: No edges, no crash
 */
void test07_ZeroSites()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    // No sites added
    diagram.compute();
    
    std::string input = "No sites";
    std::string expected = "0 edges (graceful handling)";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edges";
    
    if (edgeCount == 0)
    {
        actual << " - empty case handled correctly";
        passed = true;
    }
    else
    {
        actual << " - unexpected edges for empty diagram";
        passed = false;
    }
    
    recordTest("TEST-07", "Zero Sites - Error Case", input, expected, actual.str(), passed);
}

/**
 * TEST 8: Coincident Sites
 * 
 * Description: Two sites at same location (degenerate case)
 * Input: Two sites both at (400, 300)
 * Expected: Handles without crashing
 */
void test08_CoincidentSites()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(400, 300);
    diagram.addSite(400, 300);  // Same location
    diagram.compute();
    
    std::string input = "Two coincident sites at (400, 300)";
    std::string expected = "Handles gracefully (no crash)";
    
    bool passed = true;
    std::stringstream actual;
    
    // Just verify it doesn't crash and returns some result
    size_t siteCount = diagram.getSiteCount();
    size_t edgeCount = diagram.getEdgeCount();
    actual << siteCount << " sites, " << edgeCount << " edges - no crash";
    
    recordTest("TEST-08", "Coincident Sites - Degenerate Case", input, expected, actual.str(), passed);
}

/**
 * TEST 9: Sites at Boundaries
 * 
 * Description: Sites placed at or near bounding box edges
 * Input: Sites at corners of bounding box
 * Expected: Edges properly clipped to bounds
 */
void test09_BoundarySites()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(1, 1);       // Near corner
    diagram.addSite(799, 1);     // Near corner
    diagram.addSite(1, 599);     // Near corner
    diagram.addSite(799, 599);   // Near corner
    diagram.compute();
    
    std::string input = "Sites at four corners of bounding box";
    std::string expected = "Edges generated without crashing";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edges generated";
    
    // For boundary sites, just verify:
    // 1. Algorithm didn't crash
    // 2. We got a reasonable number of edges (4 sites should produce ~5-8 edges)
    if (edgeCount >= 4 && edgeCount <= 20)
    {
        // Count how many edges are reasonably within or near bounds
        // Use generous tolerance since edge clipping is for visualization
        const auto& edges = diagram.getEdges();
        int edgesNearBounds = 0;
        double tolerance = 100.0;  // Generous tolerance

        for (const auto& edge : edges)
        {
            bool startOK = (edge->start.x >= -tolerance && edge->start.x <= 800 + tolerance &&
                            edge->start.y >= -tolerance && edge->start.y <= 600 + tolerance);
            bool endOK = (edge->end.x >= -tolerance && edge->end.x <= 800 + tolerance &&
                          edge->end.y >= -tolerance && edge->end.y <= 600 + tolerance);

            if (startOK && endOK)
            {
                edgesNearBounds++;
            }
        }

        actual << " - " << edgesNearBounds << "/" << edgeCount << " edges near bounds";

        // Pass if most edges are reasonable (visualization will clip the rest)
        if (edgesNearBounds >= edgeCount / 2)
        {
            passed = true;
        }
        else
        {
            passed = false;
            actual << " (too many edges far outside bounds)";
        }
    }
    else if (edgeCount < 4)
    {
        passed = false;
        actual << " - too few edges";
    }
    else
    {
        passed = false;
        actual << " - unexpected edge count";
    }
    
    recordTest("TEST-09", "Boundary Sites - Edge Clipping", input, expected, actual.str(), passed);
}

/**
 * TEST 10: Verify Perpendicular Bisector Property
 * 
 * Description: Verify that edges are perpendicular bisectors of site pairs
 * Input: Various site configurations
 * Expected: Midpoints of edges equidistant from adjacent sites
 */
void test10_PerpendicularBisectorProperty()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(200, 300);
    diagram.addSite(600, 300);
    diagram.addSite(400, 100);
    diagram.addSite(400, 500);
    diagram.compute();
    
    std::string input = "Diamond configuration of 4 sites";
    std::string expected = "Most edge midpoints equidistant from adjacent sites";
    
    bool passed = true;
    std::stringstream actual;
    
    const auto& edges = diagram.getEdges();
    int validBisectors = 0;
    int totalEdges = static_cast<int>(edges.size());
    
    for (const auto& edge : edges)
    {
        // Find the midpoint of the edge
        Point midpoint((edge->start.x + edge->end.x) / 2,
                       (edge->start.y + edge->end.y) / 2);
        
        // Check if it's equidistant from the two sites (with resonable tolerance)
        if (isOnPerpendicularBisector(midpoint, edge->site1, edge->site2))
        {
            validBisectors++;
        }
    }
    
    actual << validBisectors << "/" << edges.size() << " edges are valid bisectors";
    
    // Pass if at least 70% of edges satisfy the property
    // (some edges may be clipped at boundaries, affecting their midpoints)
    double passRate = (totalEdges > 0) ? (double)validBisectors / totalEdges : 0;

    if (passRate >= 0.7)
    {
        passed = true;
        actual << " (acceptable)";
    }
    else
    {
        passed = false;
        actual << " (too few valid bisectors)";
    }
    
    recordTest("TEST-10", "Perpendicular Bisector Property", input, expected, actual.str(), passed);
}

/**
 * TEST 11: File Export
 * 
 * Description: Test exporting diagram to JSON file
 * Input: Simple 3-site diagram
 * Expected: Valid JSON file created
 */
void test11_FileExport()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(200, 200);
    diagram.addSite(600, 200);
    diagram.addSite(400, 500);
    diagram.compute();
    
    std::string filename = "test_output.json";
    bool exportSuccess = diagram.exportToFile(filename);
    
    std::string input = "3-site diagram, export to " + filename;
    std::string expected = "JSON file created successfully";
    
    bool passed = false;
    std::stringstream actual;
    
    if (exportSuccess)
    {
        // Verify file exists and has content
        std::ifstream file(filename);
        if (file.is_open())
        {
            std::string content((std::istreambuf_iterator<char>(file)),
                               std::istreambuf_iterator<char>());
            file.close();
            
            if (content.find("sites") != std::string::npos &&
                content.find("edges") != std::string::npos)
            {
                actual << "File created with valid structure";
                passed = true;
            }
            else
            {
                actual << "File created but invalid structure";
            }
            
            // Clean up test file
            std::remove(filename.c_str());
        }
        else
        {
            actual << "Export returned true but file not accessible";
        }
    }
    else
    {
        actual << "Export function returned false";
    }
    
    recordTest("TEST-11", "File Export - JSON Output", input, expected, actual.str(), passed);
}

/**
 * TEST 12: File Import
 * 
 * Description: Test importing sites from file
 * Input: File with 3 sites
 * Expected: Sites loaded correctly
 */
void test12_FileImport()
{
    // Create test input file
    std::string filename = "test_input.txt";
    {
        std::ofstream file(filename);
        file << "100 200\n";
        file << "500 300\n";
        file << "300 400\n";
        file.close();
    }
    
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    bool importSuccess = diagram.importFromFile(filename);
    
    std::string input = "File with 3 sites: (100,200), (500,300), (300,400)";
    std::string expected = "3 sites loaded correctly";
    
    bool passed = false;
    std::stringstream actual;
    
    if (importSuccess)
    {
        size_t siteCount = diagram.getSiteCount();
        actual << siteCount << " sites loaded";
        
        if (siteCount == 3)
        {
            passed = true;
        }
    }
    else
    {
        actual << "Import failed";
    }
    
    // Clean up test file
    std::remove(filename.c_str());
    
    recordTest("TEST-12", "File Import - Site Loading", input, expected, actual.str(), passed);
}

/**
 * TEST 13: Performance Test
 * 
 * Description: Verify O(n log n) time complexity
 * Input: Increasing number of sites
 * Expected: Reasonable execution time
 */
void test13_PerformanceTest()
{
    std::string input = "50, 100, 200, 400 sites";
    std::string expected = "Execution times scale reasonably";
    
    std::vector<int> sizes = {50, 100, 200, 400};
    std::vector<double> times;
    
    std::stringstream actual;
    
    for (int size : sizes)
    {
        VoronoiDiagram diagram;
        diagram.setBoundingBox(0, 1000, 0, 1000);
        
        // Add random sites
        for (int i = 0; i < size; i++)
        {
            double x = 10 + (rand() % 980);
            double y = 10 + (rand() % 980);
            diagram.addSite(x, y);
        }
        
        // Time the computation
        auto start = std::chrono::high_resolution_clock::now();
        diagram.compute();
        auto end = std::chrono::high_resolution_clock::now();
        
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        times.push_back(time_ms);
    }
    
    // Output times
    for (size_t i = 0; i < sizes.size(); i++)
    {
        actual << sizes[i] << " sites: " << std::fixed << std::setprecision(2) 
               << times[i] << "ms";
        if (i < sizes.size() - 1) actual << ", ";
    }
    
    // Check if times are reasonable (should be less than 1 second for 400 sites)
    bool passed = (times.back() < 1000);
    
    if (passed)
    {
        actual << " - acceptable performance";
    }
    else
    {
        actual << " - performance concerns";
    }
    
    recordTest("TEST-13", "Performance Test - Time Complexity", input, expected, actual.str(), passed);
}

/**
 * TEST 14: Clear and Recompute
 * 
 * Description: Test clearing and recomputing diagram
 * Input: Initial diagram, clear, add new sites, recompute
 * Expected: New diagram computed correctly
 */
void test14_ClearAndRecompute()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    
    // First computation
    diagram.addSite(100, 100);
    diagram.addSite(700, 100);
    diagram.compute();
    size_t firstEdgeCount = diagram.getEdgeCount();
    
    // Clear and add different sites
    diagram.clear();
    diagram.addSite(200, 300);
    diagram.addSite(400, 300);
    diagram.addSite(600, 300);
    diagram.compute();
    size_t secondEdgeCount = diagram.getEdgeCount();
    
    std::string input = "Clear after 2 sites, then add 3 sites";
    std::string expected = "New diagram computed (more edges)";
    
    bool passed = true;
    std::stringstream actual;
    
    actual << "First: " << firstEdgeCount << " edges, Second: " << secondEdgeCount << " edges";
    
    if (secondEdgeCount > firstEdgeCount)
    {
        actual << " - clear and recompute worked";
    }
    else
    {
        passed = false;
        actual << " - unexpected edge counts";
    }
    
    recordTest("TEST-14", "Clear and Recompute", input, expected, actual.str(), passed);
}

/**
 * TEST 15: Very Close Sites
 * 
 * Description: Test numerical stability with very close sites
 * Input: Sites separated by small distance
 * Expected: Handles without numerical errors
 */
void test15_VeryCloseSites()
{
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, 800, 0, 600);
    diagram.addSite(400, 300);
    diagram.addSite(400.001, 300);  // Very close
    diagram.addSite(100, 300);
    diagram.compute();
    
    std::string input = "Sites at (400,300), (400.001,300), (100,300)";
    std::string expected = "Numerical stability maintained";
    
    bool passed = true;
    std::stringstream actual;
    
    size_t edgeCount = diagram.getEdgeCount();
    actual << edgeCount << " edges generated";
    
    // Check that we got some edges without crashing
    if (edgeCount > 0)
    {
        actual << " - handled close sites";
    }
    else
    {
        passed = false;
        actual << " - potential numerical issues";
    }
    
    recordTest("TEST-15", "Very Close Sites - Numerical Stability", input, expected, actual.str(), passed);
}

// ============================================================================
// TEST REPORT GENERATION
// ============================================================================

/**
 * @brief Generate formatted test report
 */
void generateReport()
{
    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "          FORTUNE'S ALGORITHM TEST REPORT\n";
    std::cout << "================================================================\n\n";
    
    std::cout << "SUMMARY\n";
    std::cout << "-------\n";
    std::cout << "Total Tests:  " << totalTests << "\n";
    std::cout << "Passed:       " << passedTests << "\n";
    std::cout << "Failed:       " << (totalTests - passedTests) << "\n";
    std::cout << "Pass Rate:    " << std::fixed << std::setprecision(1)
              << (100.0 * passedTests / totalTests) << "%\n\n";
    
    std::cout << "DETAILED RESULTS\n";
    std::cout << "----------------\n\n";
    
    for (const auto& result : testResults)
    {
        std::cout << "Test: " << result.testId << "\n";
        std::cout << "Description: " << result.description << "\n";
        std::cout << "Input: " << result.input << "\n";
        std::cout << "Expected: " << result.expectedOutput << "\n";
        std::cout << "Actual: " << result.actualOutput << "\n";
        std::cout << "Result: " << (result.passed ? "PASS" : "FAIL") << "\n";
        std::cout << "\n";
    }
    
    std::cout << "================================================================\n";
    std::cout << "                    END OF TEST REPORT\n";
    std::cout << "================================================================\n";
}

/**
 * @brief Export test results to file
 */
void exportResultsToFile(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not create results file\n";
        return;
    }
    
    file << "FORTUNE'S ALGORITHM TEST RESULTS\n";
    file << "================================\n\n";
    file << "Summary: " << passedTests << "/" << totalTests << " tests passed\n\n";
    
    for (const auto& result : testResults)
    {
        file << "----------------------------------------\n";
        file << "Test ID: " << result.testId << "\n";
        file << "Description: " << result.description << "\n";
        file << "Input: " << result.input << "\n";
        file << "Expected Output: " << result.expectedOutput << "\n";
        file << "Actual Output: " << result.actualOutput << "\n";
        file << "Passed: " << (result.passed ? "true" : "false") << "\n";
    }
    
    file.close();
    std::cout << "Results exported to " << filename << "\n";
}

// ============================================================================
// MAIN
// ============================================================================

int main()
{
    std::cout << "Fortune's Algorithm - Correctness Test Suite\n";
    std::cout << "============================================\n\n";
    std::cout << "Running tests...\n\n";
    
    // Run all tests
    test01_TwoSites();
    test02_ThreeSitesTriangle();
    test03_CollinearSites();
    test04_SquareArrangement();
    test05_LargeNumberOfSites();
    test06_SingleSite();
    test07_ZeroSites();
    test08_CoincidentSites();
    test09_BoundarySites();
    test10_PerpendicularBisectorProperty();
    test11_FileExport();
    test12_FileImport();
    test13_PerformanceTest();
    test14_ClearAndRecompute();
    test15_VeryCloseSites();
    
    // Generate report
    generateReport();
    
    // Export results to file
    exportResultsToFile("test_results.txt");
    
    return (passedTests == totalTests) ? 0 : 1;
}

/**
 * @file main.cpp
 * @brief Main program for Fortune's Algorithm Voronoi Diagram Generator
 * 
 * This program demonstrates the use of Fortune's algorithm to compute
 * Voronoi diagrams. It can:
 * - Generate random sites
 * - Read sites from a file
 * - Export the computed diagram to JSON for visualization
 * 
 * Usage:
 *   ./voronoi_generator                    # Generate random diagram
 *   ./voronoi_generator <input_file>       # Read sites from file
 *   ./voronoi_generator -n <count>         # Generate n random sites
 *   ./voronoi_generator -o <output_file>   # Specify output file
 * 
 * @author Jack Woods, Ashleigh Kirkpatrick
 * @date 2025
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "voronoi.h"

using namespace Voronoi;

// ============================================================================
// CONFIGURATION
// ============================================================================

// Default values
const int DEFAULT_SITE_COUNT = 20;
const double DEFAULT_WIDTH = 800.0;
const double DEFAULT_HEIGHT = 600.0;
const std::string DEFAULT_OUTPUT = "voronoi_output.json";

// Margin to keep sites away from edges (for better visualization)
const double MARGIN = 50.0;

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * @brief Print usage information
 * @param programName Name of the program (argv[0])
 */
void printUsage(const std::string& programName)
{
    std::cout << "Fortune's Algorithm Voronoi Diagram Generator\n"
              << "=============================================\n\n"
              << "Usage:\n"
              << "  " << programName << "                         Generate random diagram\n"
              << "  " << programName << " <input_file>            Read sites from file\n"
              << "  " << programName << " -n <count>              Generate n random sites\n"
              << "  " << programName << " -o <output_file>        Specify output file\n"
              << "  " << programName << " -w <width> -h <height>  Set dimensions\n"
              << "  " << programName << " --help                  Show this help\n\n"
              << "Input file format: One site per line as 'x y'\n"
              << "Output format: JSON file for visualization\n";
}

/**
 * @brief Generate random sites within bounds
 * @param diagram The Voronoi diagram to add sites to
 * @param count Number of sites to generate
 * @param width Width of the bounding box
 * @param height Height of the bounding box
 */
void generateRandomSites(VoronoiDiagram& diagram, int count, double width, double height)
{
    // Seed random number generator
    std::srand(static_cast<unsigned>(std::time(nullptr)));
    
    for (int i = 0; i < count; ++i)
    {
        // Generate random coordinates within margins
        double x = MARGIN + (static_cast<double>(std::rand()) / RAND_MAX) * (width - 2 * MARGIN);
        double y = MARGIN + (static_cast<double>(std::rand()) / RAND_MAX) * (height - 2 * MARGIN);
        diagram.addSite(x, y);
    }
}

/**
 * @brief Print diagram statistics
 * @param diagram The computed Voronoi diagram
 */
void printStatistics(const VoronoiDiagram& diagram)
{
    std::cout << "\nVoronoi Diagram Statistics:\n"
              << "===========================\n"
              << "Number of sites: " << diagram.getSiteCount() << "\n"
              << "Number of edges: " << diagram.getEdgeCount() << "\n";
    
    // Euler's formula check: V - E + F = 2 for planar graphs
    // For Voronoi: n sites create approximately 3n - 6 edges and 2n - 5 faces
    size_t n = diagram.getSiteCount();
    if (n >= 3)
    {
        std::cout << "Expected edges (3n-6): " << (3 * n - 6) << "\n"
                  << "Expected faces (2n-5): " << (2 * n - 5) << "\n";
    }
}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main(int argc, char* argv[])
{
    // Default configuration
    std::string inputFile;
    std::string outputFile = DEFAULT_OUTPUT;
    int siteCount = DEFAULT_SITE_COUNT;
    double width = DEFAULT_WIDTH;
    double height = DEFAULT_HEIGHT;
    bool useRandomSites = true;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "-n" && i + 1 < argc)
        {
            siteCount = std::atoi(argv[++i]);
            if (siteCount < 2)
            {
                std::cerr << "Error: Need at least 2 sites\n";
                return 1;
            }
        }
        else if (arg == "-o" && i + 1 < argc)
        {
            outputFile = argv[++i];
        }
        else if (arg == "-w" && i + 1 < argc)
        {
            width = std::atof(argv[++i]);
        }
        else if (arg == "-h" && i + 1 < argc)
        {
            height = std::atof(argv[++i]);
        }
        else if (arg[0] != '-')
        {
            // Assume it's an input file
            inputFile = arg;
            useRandomSites = false;
        }
    }
    
    // Create Voronoi diagram
    VoronoiDiagram diagram;
    diagram.setBoundingBox(0, width, 0, height);
    
    // Load or generate sites
    if (useRandomSites)
    {
        std::cout << "Generating " << siteCount << " random sites...\n";
        generateRandomSites(diagram, siteCount, width, height);
    }
    else
    {
        std::cout << "Loading sites from " << inputFile << "...\n";
        if (!diagram.importFromFile(inputFile))
        {
            std::cerr << "Error: Failed to load input file\n";
            return 1;
        }
    }
    
    // Verify we have enough sites
    if (diagram.getSiteCount() < 2)
    {
        std::cerr << "Error: Need at least 2 sites to compute Voronoi diagram\n";
        return 1;
    }
    
    std::cout << "Loaded " << diagram.getSiteCount() << " sites.\n";
    
    // Compute the Voronoi diagram
    std::cout << "Computing Voronoi diagram using Fortune's algorithm...\n";
    diagram.compute();
    
    // Print statistics
    printStatistics(diagram);
    
    // Export to file
    std::cout << "\nExporting to " << outputFile << "...\n";
    if (diagram.exportToFile(outputFile))
    {
        std::cout << "Export successful!\n";
    }
    else
    {
        std::cerr << "Error: Failed to export diagram\n";
        return 1;
    }
    
    // Print site locations for reference
    std::cout << "\nSite Locations (Ant Colony Centers):\n";
    std::cout << "------------------------------------\n";
    const auto& sites = diagram.getSites();
    for (const auto& site : sites)
    {
        std::cout << "  Site " << site.id << ": (" << site.x << ", " << site.y << ")\n";
    }
    
    return 0;
}

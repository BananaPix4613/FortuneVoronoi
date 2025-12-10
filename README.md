# Fortune's Algorithm Voronoi Diagram Project

## Ant Colony Territory Visualization

**Course:** CSC 372 - Analysis of Algorithms
**Date:** 12/9/2025

---

## Project Overview

This project implements **Fortune's Algorithm** for computing Voronoi diagrams in C++, with a Python/Pygame visualization that presents the diagram as **ant colony territory boundaries**.

### What is a Voronoi Diagram?

A Voronoi diagram partitions a plane into regions based on the distance to a set of points (called "sites"). Each region contains all points closer to its site than to any other site. In our visualization:

- **Sites** = Ant colony locations (anthills)
- **Voronoi edges** = Territory boundaries where two colonies have equal influence
- **Voronoi regions** = Each colony's territory

### Fortune's Algorithm

Fortune's Algorithm is an efficient sweep line algorithm for computing Voronoi diagrams:

- **Time Complexity:** O(n log n)
- **Space Complexity:** O(n)
- **Method:** Maintains a "beach line" of parabolas as a sweep line moves across the plane

---

## Project Structure

```
fortune_voronoi/
├── include/
│   └── voronoi.h           # Header file with data structures
├── src/
│   ├── voronoi.cpp         # Fortune's algorithm implementation
│   └── main.cpp            # Main program
├── tests/
│   └── test_voronoi.cpp    # Correctness test suite
├── visualization/
│   └── ant_territory_visualizer.py  # Pygame visualization
├── VoronoiGenerator/       # Visual Studio subproject for the generator
├── TestVoronoi/            # Visual studio subproject for the test cases
├── FortuneVoronoi.sln      # Visual Studio project solution file
├── LICENSE                 # License details
└── README.md               # This file
```

---

## Building the Project

### Prerequisites

**Windows:**
- **Visual Studio C++ Compiler:** Visual Studio with C++17 support
- **Python 3.x** with Pygame (`pip install pygame`)

### Build Commands

**Visual Studio**

Open project through the .sln solution file and click build from the dropdown menu from the top of the application.

## Running the Project

### Step 1: Generate a Voronoi Diagram

**Windows:**
```cmd
:: Run directly
VoronoiGenerator.exe

:: Or with parameters
VoronoiGenerator.exe -n 30
```

### Step 2: Visualize as Ant Territories

**Windows:**
```cmd
python visualization\ant_territory_visualizer.py
```

---

## Visualization Controls

|  Key  |          Action          |
|-------|--------------------------|
| SPACE | Skip sweep animation     |
|   P   | Pause/Resume animation   |
|   A   | Toggle ant display       |
|   B   | Toggle boundary display  |
|   R   | Reset ants to colonies   |
|   C   | Sweep animation replay   |
|   Q   | Quit visualization       |
| Click | Select a colony for info |

---

## Running Tests

The test suite verifies the correctness of the Fortune's algorithm implementation.

```cmd
:: Run directly
TestVoronoi.exe
```

### Test Categories

1. **Normal Usage** - Basic diagram generation
2. **Edge Cases** - Two sites, collinear sites, etc.
3. **Error Handling** - Zero sites, single site
4. **Numerical Stability** - Very close sites
5. **Performance** - Large number of sites

---

## Input File Format

Sites can be loaded from a text file with one point per line:

```
# Comments start with #
100 200
500 300
300 450
# x y coordinates separated by whitespace
```

---

## Output File Format (JSON)

The generated diagram is saved in JSON format:

```json
{
  "bounds": {
    "minX": 0, "maxX": 800,
    "minY": 0, "maxY": 600
  },
  "sites": [
    {"x": 100, "y": 200, "id": 0},
    {"x": 500, "y": 300, "id": 1}
  ],
  "edges": [
    {
      "start": {"x": 300, "y": 0},
      "end": {"x": 300, "y": 600},
      "site1_id": 0,
      "site2_id": 1
    }
  ]
}
```

---

## Algorithm Details

### Key Data Structures

1. **Point** - 2D coordinate with ID
2. **Edge** - Line segment between two vertices (bisector of two sites)
3. **Arc** - Parabolic arc in the beach line
4. **Event** - Site event or circle event

### Algorithm Steps

1. **Initialize:** Add all sites to priority queue as site events
2. **Process Events:** While queue not empty:
   - **Site Events:** Insert new arc, create new edges
   - **Circle Event:** Remove arc, complete vertex
3. **Finish:** Clip infinite edges to bounding box

### Complexity Analysis

|  Operation   | Complexity |
|--------------|------------|
|   Overall    | O(n log n) |
|  Site Event  |  O(log n)  |
| Circle Event |  O(log n)  |
|    Space     |    O(n)    |

---

## Test Plan Summary

See `tests/test_voronoi.cpp` for complete test implementation.

| Test ID |      Description       |   Category   |
|---------|------------------------|--------------|
| TEST-01 | Two sites (minimum)    | Normal       |
| TEST-02 | Three sites (triangle) | Normal       |
| TEST-03 | Collinear sites        | Edge Case    |
| TEST-04 | Square arrangement     | Normal       |
| TEST-05 | Large site count       | Stress       |
| TEST-06 | Single site            | Error        |
| TEST-07 | Zero sites             | Error        |
| TEST-08 | Coincident sites       | Error        |
| TEST-09 | Boundary sites         | Edge Case    |
| TEST-10 | Bisector property      | Verification |
| TEST-11 | File export            | I/O          |
| TEST-12 | File import            | I/O          |
| TEST-13 | Performance            | Performance  |
| TEST-14 | Clear/recompute        | Normal       |
| TEST-15 | Close sites            | Numerical    |

---

## Acknowledgments

- Fortune's algorithm: Steven Fortune (1986)
- Pygame visualization framework
- CSC 372 course materials

---

## Authors

Jack Woods, Ashleigh Kirkpatrick - 2025
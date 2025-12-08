#!/usr/bin/env python3
"""
ant_territory_visualizer.py
============================
Pygame visualization of Voronoi diagrams as ant colony territory boundaries.

This script reads a JSON file produced by the C++ Fortune's algorithm implementation
and displays it as a visual representation of ant colony territories. Each site
represents an ant colony, and the Voronoi edges represent the territorial boundaries
where two colonies' influence meets.

Features:
- Animated sweep line showing Fortune's algorithm in action
- Territory-bounded ants that stay within their Voronoi cells
- Colorful territory regions for each ant colony
- Interactive controls

Author: Jack Woods
Date: 2025

Requirements:
- Python 3.x
- Pygame (pip install pygame)
"""

import pygame
import json
import sys
import random
import math
from typing import List, Tuple, Optional
from dataclasses import dataclass
from enum import Enum, auto

# ============================================================================
# CONFIGURATION
# ============================================================================

# Window settings
WINDOW_WIDTH = 800
WINDOW_HEIGHT = 600
FPS = 60

# Colors (RGB tuples)
BACKGROUND_COLOR = (255, 255, 255)     # White background for sweep animation
SIMULATION_BG_COLOR = (34, 139, 34)    # Forest green for ant simulation
BOUNDARY_COLOR = (139, 69, 19)         # Saddle brown for boundaries
SITE_COLOR = (0, 0, 0)                 # Black for sites
TEXT_COLOR = (255, 255, 255)
SWEEP_LINE_COLOR = (255, 0, 0)         # Red sweep line
BEACH_LINE_COLOR = (0, 100, 255)       # Blue for beach line
EDGE_COLOR = (100, 100, 100)           # Gray for edges during animation

# Ant territory colors (pastel palette for visibility)
TERRITORY_COLORS = [
    (255, 182, 193),   # Light pink
    (255, 218, 185),   # Peach
    (255, 255, 186),   # Light yellow
    (186, 255, 201),   # Mint green
    (186, 225, 255),   # Light blue
    (221, 186, 255),   # Lavender
    (255, 186, 255),   # Light magenta
    (186, 255, 255),   # Light cyan
    (255, 228, 196),   # Bisque
    (240, 230, 140),   # Khaki
]

# Animation settings
NUM_ANTS_PER_COLONY = 5
ANT_SPEED = 1.5
ANT_SIZE = 4
SWEEP_LINE_SPEED = 2  # Pixels per frame during algorithm animation

# UI settings
INFO_PANEL_HEIGHT = 80
TITLE = "Ant Colony Territory Map - Fortune's Algorithm Visualization"


# ============================================================================
# VISUALIZATION STATES
# ============================================================================

class VisualizationState(Enum):
    """States of the visualization"""
    SWEEP_ANIMATION = auto()    # Showing Fortune's algorithm
    ANT_SIMULATION = auto()     # Normal ant territory view


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class Point:
    """Represents a 2D point"""
    x: float
    y: float
    
    def to_tuple(self) -> Tuple[int, int]:
        """Convert to integer tuple for Pygame drawing"""
        return (int(self.x), int(self.y))
    
    def distance_to(self, other: 'Point') -> float:
        """Calculate Euclidean distance to another point"""
        return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)


@dataclass
class Site:
    """Represents an ant colony location"""
    position: Point
    id: int
    color: Tuple[int, int, int]
    name: str


@dataclass
class Edge:
    """Represents a territorial boundary (Voronoi edge)"""
    start: Point
    end: Point
    site1_id: int
    site2_id: int


@dataclass
class Ant:
    """Represents an animated ant in the visualization"""
    position: Point
    velocity: Point
    colony_id: int
    color: Tuple[int, int, int]

    def update(self, sites: List[Site], bounds: Tuple[float, float, float, float]):
        """
        Update ant position, keeping it within its home territory.
        
        Args:
            sites: List of all colony sites
            edges: List of Voronoi edges (territory boundaries)
            bounds: (min_x, max_x, min_y, max_y) bounding box
        """
        # Find home colony position
        home_colony = None
        for site in sites:
            if site.id == self.colony_id:
                home_colony = site
                break

        if home_colony is None:
            return

        # Add random wandering
        self.velocity.x += random.uniform(-0.3, 0.3)
        self.velocity.y += random.uniform(-0.3, 0.3)
        
        # Limit speed
        speed = math.sqrt(self.velocity.x**2 + self.velocity.y**2)
        if speed > ANT_SPEED:
            self.velocity.x = (self.velocity.x / speed) * ANT_SPEED
            self.velocity.y = (self.velocity.y / speed) * ANT_SPEED
        
        # Calculate next position
        next_x = self.position.x + self.velocity.x
        next_y = self.position.y + self.velocity.y
        next_pos = Point(next_x, next_y)

        # # Check if next position would be in a different territory
        # next_owner = self._get_closest_site_id(next_pos, sites)

        # if next_owner != self.colony_id:
        #     # Would leave territory - bounce back toward home
        #     dx = home_colony.position.x - self.position.x
        #     dy = home_colony.position.y - self.position.y
        #     dist = math.sqrt(dx**2 + dy**2)

        #     if dist > 0:
        #         # Reflect velocity and add push toward home
        #         self.velocity.x = (dx / dist) * ANT_SPEED * 0.8
        #         self.velocity.y = (dy / dist) * ANT_SPEED * 0.8
            
        #     # Don't update position - stay in current spot
        #     return

        # Safe to move - update position
        self.position.x = next_x
        self.position.y = next_y
        
        # Gentle attraction to colony center when far away
        dx = home_colony.position.x - self.position.x
        dy = home_colony.position.y - self.position.y
        dist = math.sqrt(dx**2 + dy**2)

        if dist > 80:
            self.velocity.x += (dx / dist) * 0.1
            self.velocity.y += (dy / dist) * 0.1
        
        # Bounce off world boundaries
        if self.position.x < bounds[0] + 5:
            self.position.x = bounds[0] + 5
            self.velocity.x = abs(self.velocity.x)
        elif self.position.x > bounds[1] - 5:
            self.position.x = bounds[1] - 5
            self.velocity.x = -abs(self.velocity.x)
            
        if self.position.y < bounds[2] + 5:
            self.position.y = bounds[2] + 5
            self.velocity.y = abs(self.velocity.y)
        elif self.position.y > bounds[3] - 5:
            self.position.y = bounds[3] - 5
            self.velocity.y = -abs(self.velocity.y)
    
    def _get_closest_site_id(self, pos: Point, sites: List[Site]) -> int:
        """Find the ID of the closest site to a position"""
        min_dist = float('inf')
        closest_id = self.colony_id
        
        for site in sites:
            dist = pos.distance_to(site.position)
            if dist < min_dist:
                min_dist = dist
                closest_id = site.id
        
        return closest_id


# ============================================================================
# VORONOI DATA LOADER
# ============================================================================

class VoronoiData:
    """
    Loads and stores Voronoi diagram data from JSON file.
    
    The JSON file should have the structure:
    {
        "bounds": {"minX": ..., "maxX": ..., "minY": ..., "maxY": ...},
        "sites": [{"x": ..., "y": ..., "id": ...}, ...],
        "edges": [{"start": {"x": ..., "y": ...}, "end": {"x": ..., "y": ...}, 
                   "site1_id": ..., "site2_id": ...}, ...]
    }
    """
    
    def __init__(self):
        self.sites: List[Site] = []
        self.edges: List[Edge] = []
        self.bounds: Tuple[float, float, float, float] = (0, 800, 0, 600)
    
    def load_from_file(self, filename: str) -> bool:
        """
        Load Voronoi data from a JSON file.
        
        Args:
            filename: Path to the JSON file
            
        Returns:
            True if loading was successful, False otherwise
        """
        try:
            with open(filename, 'r') as f:
                data = json.load(f)
            
            # Load bounds
            if 'bounds' in data:
                b = data['bounds']
                self.bounds = (b['minX'], b['maxX'], b['minY'], b['maxY'])
            
            # Load sites
            self.sites = []
            for i, site_data in enumerate(data.get('sites', [])):
                color = TERRITORY_COLORS[i % len(TERRITORY_COLORS)]
                site = Site(
                    position=Point(site_data['x'], site_data['y']),
                    id=site_data.get('id', i),
                    color=color,
                    name=f"Colony {site_data.get('id', i)}"
                )
                self.sites.append(site)
            
            # Load edges
            self.edges = []
            for edge_data in data.get('edges', []):
                start = Point(edge_data['start']['x'], edge_data['start']['y'])
                end = Point(edge_data['end']['x'], edge_data['end']['y'])
                
                # Filter out invalid edges (starting at origin or with invalid coordinates)
                if (abs(start.x) < 1 and abs(start.y) < 1):
                    continue  # Skip edges starting at (0,0) - uninitialized
                if start.x < -1000 or start.y < -1000 or end.x < -1000 or end.y < -1000:
                    continue  # Skip edges with extreme negative coords
                if start.x > 10000 or start.y > 10000 or end.x > 10000 or end.y > 10000:
                    continue  # Skip edges with extreme positive coords
                    
                edge = Edge(
                    start=start,
                    end=end,
                    site1_id=edge_data.get('site1_id', -1),
                    site2_id=edge_data.get('site2_id', -1)
                )
                self.edges.append(edge)
            
            print(f"Loaded {len(self.sites)} colonies and {len(self.edges)} boundaries")
            return True
            
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found")
            return False
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in '{filename}': {e}")
            return False
        except Exception as e:
            print(f"Error loading data: {e}")
            return False


# ============================================================================
# SWEEP LINE ANIMATOR
# ============================================================================

class SweepLineAnimator:
    """
    Animates Fortune's algorithm sweep line process.

    Shows:
    - Sweep line moving from top to bottom
    - Sites being "discovered" as sweep line passes
    - Edges being drawn progressively
    - Beach line parabolas (simplified visualization)
    """

    def __init__(self, data: VoronoiData, display_height: int):
        self.data = data
        self.display_height = display_height
        self.width = int(data.bounds[1] - data.bounds[0])

        # Sort sites by y-coordinate (top to bottom for animation)
        self.sorted_sites = sorted(data.sites, key=lambda s: s.position.y)

        # Animation state
        self.sweep_y = 0  # Current sweep line position
        self.discovered_sites: List[Site] = []
        self.animation_complete = False

        # Track which edges are "active" (being traced)
        # An edge is active when both its adjacent sites are discovered
        self.active_edges: List[Edge] = []
    
    def _get_parabola_y(self, focus_x: float, focus_y: float, directrix_y: float, x: float) -> Optional[float]:
        """
        Calculate y-coordinate on parabola at given x.
        
        For a parabola with focus (fx, fy) and directrix y = dy:
        The parabola is the set of points equidistant from focus and directrix.
        
        When directrix_y > focus_y (sweep line below focus), parabola opens upward.
        """
        # Distance from focus to directrix
        p = directrix_y - focus_y
        
        if abs(p) < 0.001:
            # Focus is on the directrix - degenerate case
            return None
        
        # Parabola equation: (x - fx)^2 = 2p(y - (fy + dy)/2)
        # Solving for y: y = (x - fx)^2 / (2p) + (fy + dy) / 2
        dx = x - focus_x
        y = -(dx * dx) / (2 * p) + (focus_y + directrix_y) / 2
        
        return y

    def _compute_beach_line(self) -> List[Tuple[int, int]]:
        """
        Compute the beach line as the lower envelope of all parabolas.
        
        For each x, find the MAXIMUM y value among all parabolas from discovered sites.
        This forms the beach line - the boundary closest to the sweep line.
        """
        if not self.discovered_sites or self.sweep_y <= 0:
            return []
        
        points = []
        
        for x in range(0, self.width, 3):
            max_y = -float('inf')
            
            for site in self.discovered_sites:
                # Only consider sites above the sweep line
                if site.position.y >= self.sweep_y:
                    continue
                
                y = self._get_parabola_y(
                    site.position.x, site.position.y,
                    self.sweep_y, x
                )
                
                if y is not None and y > max_y and y < self.sweep_y:
                    max_y = y
            
            if max_y > -float('inf') and max_y >= 0:
                points.append((x, int(max_y)))
        
        return points
    
    def _get_edge_visibility(self, edge: Edge) -> Tuple[bool, Optional[Point], Optional[Point]]:
        """
        Determine if an edge should be visible and how much of it.
        
        An edge becomes visible when both adjacent sites have been discovered.
        The visible portion extends from the start point to where it intersects
        with the current sweep line (or to its end point if fully visible).
        """
        # Find the sites for this edge
        site1 = None
        site2 = None
        for site in self.data.sites:
            if site.id == edge.site1_id:
                site1 = site
            if site.id == edge.site2_id:
                site2 = site
        
        if site1 is None or site2 is None:
            return False, None, None
        
        # Edge is only visible if both sites have been discovered
        if site1 not in self.discovered_sites or site2 not in self.discovered_sites:
            return False, None, None
        
        # Calculate which part of the edge is visible
        start = edge.start
        end = edge.end
        
        # Determine the "top" point (smaller y) and "bottom" point (larger y)
        if start.y <= end.y:
            top_point = start
            bottom_point = end
        else:
            top_point = end
            bottom_point = start
        
        # If the entire edge is above sweep line, show all of it
        if bottom_point.y <= self.sweep_y:
            return True, edge.start, edge.end
        
        # If edge starts below sweep line, don't show yet
        if top_point.y > self.sweep_y:
            return False, None, None
        
        # Partial visibility - interpolate to sweep line
        if abs(bottom_point.y - top_point.y) < 0.001:
            return True, edge.start, edge.end
        
        t = (self.sweep_y - top_point.y) / (bottom_point.y - top_point.y)
        t = max(0, min(1, t))
        
        clip_x = top_point.x + t * (bottom_point.x - top_point.x)
        clip_y = self.sweep_y
        
        # Return with original orientation
        if start.y <= end.y:
            return True, start, Point(clip_x, clip_y)
        else:
            return True, Point(clip_x, clip_y), end
    
    def update(self) -> bool:
        """
        Update the sweep line animation.

        Returns:
            True if animation is still running, False if complete
        """
        if self.animation_complete:
            return False
        
        # Move sweep line down
        self.sweep_y += SWEEP_LINE_SPEED

        # Check for newly discovered sites
        for site in self.sorted_sites:
            if site not in self.discovered_sites and site.position.y <= self.sweep_y:
                self.discovered_sites.append(site)
        
        # Check if animation is complete
        if self.sweep_y >= self.display_height:
            self.animation_complete = True
            return False
        
        return True
    
    def draw(self, screen: pygame.Surface, font: pygame.font.Font):
        """Draw the current state of the sweep line animation"""
        width = screen.get_width()
        height = screen.get_height() - INFO_PANEL_HEIGHT

        # Draw background
        screen.fill(BACKGROUND_COLOR)

        # Draw all edges (partially visible based on sweep line)
        for edge in self.data.edges:
            visible, start, end = self._get_edge_visibility(edge)
            if visible and start and end:
                # Clip to screen bounds for drawing
                start_x = max(0, min(width, int(start.x)))
                start_y = max(0, min(height, int(start.y)))
                end_x = max(0, min(width, int(end.x)))
                end_y = max(0, min(height, int(end.y)))
                
                pygame.draw.line(screen, EDGE_COLOR,
                                (start_x, start_y), (end_x, end_y), 2)

        # Draw the beach line (envelope of parabolas)
        beach_points = self._compute_beach_line()
        if len(beach_points) > 1:
            pygame.draw.lines(screen, BEACH_LINE_COLOR, False, beach_points, 2)
        
        # Draw all sites (undiscovered as hollow, discovered as filled)
        for site in self.data.sites:
            pos = site.position.to_tuple()
            if site in self.discovered_sites:
                pygame.draw.circle(screen, SITE_COLOR, pos, 6)
            else:
                pygame.draw.circle(screen, SITE_COLOR, pos, 6, 2)
        
        # Draw sweep line
        pygame.draw.line(screen, SWEEP_LINE_COLOR,
                        (0, int(self.sweep_y)),
                        (width, int(self.sweep_y)), 2)


# ============================================================================
# MAIN VISUALIZATION CLASS
# ============================================================================

class AntTerritoryVisualizer:
    """
    Main visualization class using Pygame.
    
    Displays the Voronoi diagram as ant colony territories with:
    - Colored territory regions
    - Animated ants
    - Territory boundary lines
    - Information panel
    """
    
    def __init__(self, voronoi_data: VoronoiData):
        """
        Initialize the visualizer.
        
        Args:
            voronoi_data: Loaded Voronoi diagram data
        """
        self.data = voronoi_data
        self.ants: List[Ant] = []
        self.paused = False
        self.selected_colony: Optional[int] = None
        self.show_ants = True
        self.show_boundaries = True
        
        # Initialize Pygame
        pygame.init()
        pygame.display.set_caption(TITLE)
        
        # Calculate window size based on data bounds
        self.width = int(voronoi_data.bounds[1] - voronoi_data.bounds[0])
        self.height = int(voronoi_data.bounds[3] - voronoi_data.bounds[2]) + INFO_PANEL_HEIGHT
        self.screen = pygame.display.set_mode((self.width, self.height))
        
        # Create clock for consistent frame rate
        self.clock = pygame.time.Clock()
        
        # Font for text display
        self.font = pygame.font.Font(None, 24)
        self.title_font = pygame.font.Font(None, 32)
        
        # Start with sweep line animation
        self.state = VisualizationState.SWEEP_ANIMATION
        self.sweep_animator = SweepLineAnimator(voronoi_data, self.height - INFO_PANEL_HEIGHT)

        # Territory surface (pre-rendered for performance)
        self.territory_surface = None
    
    def _create_ants(self):
        """Create animated ants for each colony"""
        self.ants = []
        for site in self.data.sites:
            # Create darker version of territory color for ants
            ant_color = tuple(max(0, c - 100) for c in site.color)
            
            for _ in range(NUM_ANTS_PER_COLONY):
                # Random starting position near colony center
                offset_x = random.uniform(-50, 50)
                offset_y = random.uniform(-50, 50)
                
                ant = Ant(
                    position=Point(site.position.x + offset_x, site.position.y + offset_y),
                    velocity=Point(random.uniform(-1, 1), random.uniform(-1, 1)),
                    colony_id=site.id,
                    color=ant_color
                )
                self.ants.append(ant)
    
    def _render_territory_surface(self):
        """Pre-render territory colors for performance"""
        self.territory_surface = pygame.Surface((self.width, self.height - INFO_PANEL_HEIGHT))
        self.territory_surface.fill(BACKGROUND_COLOR)

        step = 8
        for x in range(0, self.width, step):
            for y in range(0, self.height - INFO_PANEL_HEIGHT, step):
                closest = self._get_closest_site(x, y)
                if closest:
                    rect = pygame.Rect(x, y, step, step)
                    pygame.draw.rect(self.territory_surface, closest.color, rect)

    def _get_closest_site(self, x: float, y: float) -> Optional[Site]:
        """Find the site (colony) closest to a given point"""
        closest = None
        min_dist = float('inf')
        
        for site in self.data.sites:
            dist = site.position.distance_to(Point(x, y))
            if dist < min_dist:
                min_dist = dist
                closest = site
        
        return closest
    
    def _draw_territories(self):
        """Draw pre-rendered territory surface"""
        if self.territory_surface is None:
            self._render_territory_surface()
        self.screen.blit(self.territory_surface, (0, 0))
    
    def _draw_boundaries(self):
        """Draw territory boundary lines (Voronoi edges)"""
        if not self.show_boundaries:
            return
            
        for edge in self.data.edges:
            # Clip to visible area
            start_x = max(-50, min(self.width + 50, edge.start.x))
            start_y = max(-50, min(self.height + 50, edge.start.y))
            end_x = max(-50, min(self.width + 50, edge.end.x))
            end_y = max(-50, min(self.height + 50, edge.end.y))
            
            # Skip if completely outside
            if (start_x == end_x and (start_x <= -50 or start_x >= self.width + 50)):
                continue
            if (start_y == end_y and (start_y <= -50 or start_y >= self.height + 50)):
                continue
            
            start = (int(start_x), int(start_y))
            end = (int(end_x), int(end_y))
            
            # Main line
            pygame.draw.line(self.screen, BOUNDARY_COLOR, start, end, 3)
    
    def _draw_colonies(self):
        """Draw colony center markers"""
        for site in self.data.sites:
            pos = site.position.to_tuple()
            
            # Draw colony marker (anthill representation)
            pygame.draw.circle(self.screen, (139, 69, 19), pos, 15)  # Brown outer
            pygame.draw.circle(self.screen, (101, 67, 33), pos, 12)  # Darker middle
            pygame.draw.circle(self.screen, (60, 40, 20), pos, 8)    # Dark center
            
            # Draw colony ID
            text = self.font.render(str(site.id), True, SITE_COLOR)
            text_rect = text.get_rect(center=(pos[0], pos[1] - 25))
            self.screen.blit(text, text_rect)
            
            # Highlight selected colony
            if self.selected_colony == site.id:
                pygame.draw.circle(self.screen, HIGHLIGHT_COLOR, pos, 20, 3)
    
    def _draw_ants(self):
        """Draw animated ants"""
        if not self.show_ants:
            return
            
        for ant in self.ants:
            pos = ant.position.to_tuple()
            
            # Draw ant body
            pygame.draw.circle(self.screen, ant.color, pos, ANT_SIZE)
            
            # Draw ant direction indicator
            dir_x = int(pos[0] + ant.velocity.x * 3)
            dir_y = int(pos[1] + ant.velocity.y * 3)
            pygame.draw.line(self.screen, ant.color, pos, (dir_x, dir_y), 2)
    
    def _draw_info_panel(self):
        """Draw the information panel at the bottom"""
        # Panel background
        panel_rect = pygame.Rect(0, self.height - INFO_PANEL_HEIGHT, self.width, INFO_PANEL_HEIGHT)
        pygame.draw.rect(self.screen, (40, 40, 40), panel_rect)
        pygame.draw.line(self.screen, (100, 100, 100), 
                        (0, self.height - INFO_PANEL_HEIGHT), 
                        (self.width, self.height - INFO_PANEL_HEIGHT), 2)
        
        if self.state == VisualizationState.SWEEP_ANIMATION:
            # Sweep animation info
            title = self.title_font.render("Fortune's Algorithm - Sweep Line Animation", True, TEXT_COLOR)
            self.screen.blit(title, (10, self.height - INFO_PANEL_HEIGHT + 10))

            sites_found = len(self.sweep_animator.discovered_sites)
            total_sites = len(self.data.sites)
            progress = int((self.sweep_animator.sweep_y / (self.height - INFO_PANEL_HEIGHT)) * 100)

            stats_text = f"Sites discovered: {sites_found}/{total_sites} | Progress: {progress}%"
            stats = self.font.render(stats_text, True, TEXT_COLOR)
            self.screen.blit(stats, (10, self.height - INFO_PANEL_HEIGHT + 40))

            controls_text = "[SPACE] Skip animation | [Q]uit"
            controls = self.font.render(controls_text, True, (180, 180, 180))
            self.screen.blit(controls, (10, self.height - INFO_PANEL_HEIGHT + 60))
        else:
            # Normal simulation info
            title = self.title_font.render("Ant Colony Territory Map", True, TEXT_COLOR)
            self.screen.blit(title, (10, self.height - INFO_PANEL_HEIGHT + 10))
            
            stats_text = f"Colonies: {len(self.data.sites)} | Boundaries: {len(self.data.edges)} | Ants: {len(self.ants)}"
            stats = self.font.render(stats_text, True, TEXT_COLOR)
            self.screen.blit(stats, (10, self.height - INFO_PANEL_HEIGHT + 40))
            
            controls_text = "[P]ause | [A]nts | [B]oundaries | [R]eset | [S]weep replay | [Q]uit"
            controls = self.font.render(controls_text, True, (180, 180, 180))
            self.screen.blit(controls, (10, self.height - INFO_PANEL_HEIGHT + 60))
            
            # Pause indicator
            if self.paused:
                pause_text = self.title_font.render("PAUSED", True, (255, 100, 100))
                self.screen.blit(pause_text, (self.width - 100, self.height - INFO_PANEL_HEIGHT + 10))
    
    def _update_ants(self):
        """Update ant positions"""
        if self.paused:
            return
        
        for ant in self.ants:
            ant.update(self.data.sites, self.data.bounds)
    
    def _handle_click(self, pos: Tuple[int, int]):
        """Handle mouse click - select nearest colony"""
        if pos[1] < self.height - INFO_PANEL_HEIGHT:
            closest = self._get_closest_site(pos[0], pos[1])
            if closest:
                self.selected_colony = closest.id
    
    def _transition_to_simulation(self):
        """Transition from sweep animation to ant simulation"""
        self.state = VisualizationState.ANT_SIMULATION
        self._render_territory_surface()
        self._create_ants()

    def _restart_sweep_animation(self):
        """Restart the sweep line animation"""
        self.state = VisualizationState.SWEEP_ANIMATION
        self.sweep_animator = SweepLineAnimator(self.data, self.height - INFO_PANEL_HEIGHT)
        self.ants = []

    def run(self):
        """Main visualization loop"""
        running = True
        
        while running:
            # Handle events
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
                    
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_q or event.key == pygame.K_ESCAPE:
                        running = False
                    elif event.key == pygame.K_SPACE:
                        if self.state == VisualizationState.SWEEP_ANIMATION:
                            self._transition_to_simulation()
                    elif event.key == pygame.K_p:
                        if self.state == VisualizationState.ANT_SIMULATION:
                            self.paused = not self.paused
                    elif event.key == pygame.K_a:
                        self.show_ants = not self.show_ants
                    elif event.key == pygame.K_b:
                        self.show_boundaries = not self.show_boundaries
                    elif event.key == pygame.K_r:
                        if self.state == VisualizationState.ANT_SIMULATION:
                            self._create_ants()
                    elif event.key == pygame.K_s:
                        self._restart_sweep_animation()
                        
                elif event.type == pygame.MOUSEBUTTONDOWN:
                    if event.button == 1:
                        self._handle_click(event.pos)
            
            # Update and draw based on state
            if self.state == VisualizationState.SWEEP_ANIMATION:
                # Update sweep animation
                if not self.sweep_animator.update():
                    self._transition_to_simulation()

                # Draw sweep animation
                self.sweep_animator.draw(self.screen, self.font)

            else:  # ANT_SIMULATION
                # Update ants
                self._update_ants()
                
                # Draw
                self._draw_territories()
                self._draw_boundaries()
                self._draw_colonies()
                self._draw_ants()
            
            # Always draw info panel
            self._draw_info_panel()
            
            # Update display
            pygame.display.flip()
            
            # Maintain frame rate
            self.clock.tick(FPS)
        
        pygame.quit()


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """
    Main entry point for the visualization.
    
    Usage:
        python ant_territory_visualizer.py [voronoi_output.json]
    """
    # Get input filename from command line or use default
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "voronoi_output.json"
    
    print("=" * 60)
    print("  Ant Colony Territory Visualizer")
    print("  Fortune's Algorithm Voronoi Diagram")
    print("=" * 60)
    print()
    
    # Load Voronoi data
    print(f"Loading data from: {input_file}")
    data = VoronoiData()
    if not data.load_from_file(input_file):
        print("\nPlease run the C++ Voronoi generator first:")
        print("  ./voronoi_generator")
        print("Then run this visualizer with the output file.")
        sys.exit(1)
    
    print()
    print("Controls:")
    print("  SPACE - Skip sweep animation")
    print("  P - Pause/Resume ant simulation")
    print("  A - Toggle ant display")
    print("  B - Toggle boundary display")
    print("  R - Reset ants")
    print("  S - Replay sweep animation")
    print("  Q - Quit")
    print("  Click - Select colony")
    print()
    print("Starting visualization...")
    print()
    
    # Create and run visualizer
    visualizer = AntTerritoryVisualizer(data)
    visualizer.run()
    
    print("Visualization complete!")


if __name__ == "__main__":
    main()

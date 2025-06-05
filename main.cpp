#include <tchar.h>
#include <windows.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

// Forward declaration of the menu creation function
void menu(HWND hwnd);

//-------------------------------------------------------------------------------------------------
// Data Structures and Global Variables
//-------------------------------------------------------------------------------------------------

// Represents a drawn shape and its properties
struct Shape {
    std::string type;   // Shape type ("DDALine", "CirclePolar", etc.)
    int x1, y1;         // Primary coordinates (e.g., starting point)
    int x2, y2;         // Secondary coordinates (e.g., ending point) or unused
    int a, b;           // Additional parameters (e.g., radius, ellipse axes)
    COLORREF color;     // Drawing color
};

// Global vector to store all drawn shapes
static std::vector<Shape> shapes;

// Coordinates and parameters used during mouse events
static int xc, yc;                // First click (center or start)
static int xe, ye;                // Second click (radius or end)
static int xe2, ye2;              // Third click (for polar ellipse)
static int R, R2;                 // Radii or auxiliary parameters
static int windowX, windowY, windowR;  // Clipping window center and radius

// Handle for the main menu
static HMENU hMenu;

//-------------------------------------------------------------------------------------------------
// Utility Functions (Bezier, LERP, etc.)
//-------------------------------------------------------------------------------------------------

double LERP(double start, double end, double t)
{
    return ((1 - t) * start) + (t * end);
}

double quadraticBezier(double x, double y, double z, double t)
{
    double a = LERP(x, y, t);
    double b = LERP(y, z, t);
    return LERP(a, b, t);
}

double cubicBezier(double w, double x, double y, double z, double t)
{
    double a = LERP(w, x, t);
    double b = LERP(x, y, t);
    double c = LERP(y, z, t);
    double d = LERP(a, b, t);
    double e = LERP(b, c, t);
    return LERP(d, e, t);
}

//-------------------------------------------------------------------------------------------------
// File I/O: Save and Load Shapes
//-------------------------------------------------------------------------------------------------

// Saves all shapes to a plain text file (one shape per line)
void PerformSave(const char* filename) {
    ofstream out(filename);
    if (!out.is_open()) return;

    for (const auto& s : shapes) {
        // Format: type x1 y1 x2 y2 a b color
        out << s.type << " "
            << s.x1 << " " << s.y1 << " "
            << s.x2 << " " << s.y2 << " "
            << s.a  << " " << s.b  << " "
            << s.color << "\n";
    }

    out.close();
}

// Loads shapes from a plain text file and repopulates the global vector
void PerformLoad(const char* filename) {
    ifstream in(filename);
    if (!in.is_open()) return;

    shapes.clear();  // Clear existing shapes

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);

        Shape s;
        iss >> s.type
            >> s.x1 >> s.y1
            >> s.x2 >> s.y2
            >> s.a  >> s.b
            >> s.color;

        shapes.push_back(s);
    }

    in.close();
    // Invalidate handled in calling code to trigger WM_PAINT
}

//-------------------------------------------------------------------------------------------------
// Line-Drawing Algorithms
//-------------------------------------------------------------------------------------------------

// Bresenham’s Line Algorithm
void lineBresenham(HDC hdc, int x0, int y0, int x1, int y1, COLORREF color)
{
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);
    int x = x0, y = y0;
    int sx = (x0 < x1 ? 1 : -1);
    int sy = (y0 < y1 ? 1 : -1);
    int err = dx - dy;

    while (true) {
        SetPixel(hdc, x, y, color);
        if (x == x1 && y == y1) break;
        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (e2 < dx) {
            err += dx;
            y += sy;
        }
    }
}

// DDA Line Algorithm
void DDALine(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    float Xinc = dx / (float)steps;
    float Yinc = dy / (float)steps;
    float X = static_cast<float>(x1);
    float Y = static_cast<float>(y1);

    for (int i = 0; i <= steps; i++) {
        SetPixel(hdc, static_cast<int>(round(X)), static_cast<int>(round(Y)), color);
        X += Xinc;
        Y += Yinc;
    }
}

// Parametric Line Algorithm (with optional circular clipping)
void parametricLine(HDC hdc, int x1, int y1, int x2, int y2, int cx, int cy, int radius, COLORREF color)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    if (steps == 0) {
        if (!radius || ((x1 - cx)*(x1 - cx) + (y1 - cy)*(y1 - cy) <= radius*radius)) {
            SetPixel(hdc, x1, y1, color);
        }
        return;
    }

    double tStep = 1.0 / steps;
    for (double t = 0.0; t <= 1.0; t += tStep) {
        int x = static_cast<int>(round(LERP(x1, x2, t)));
        int y = static_cast<int>(round(LERP(y1, y2, t)));
        if (!radius || ((x - cx)*(x - cx) + (y - cy)*(y - cy) <= radius*radius)) {
            SetPixel(hdc, x, y, color);
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Circle-Drawing Algorithms
//-------------------------------------------------------------------------------------------------

// Draws 8 symmetric points around the center
void draw8Points(HDC hdc, int cx, int cy, int x, int y, COLORREF color)
{
    SetPixel(hdc, cx + x, cy + y, color);
    SetPixel(hdc, cx - x, cy + y, color);
    SetPixel(hdc, cx - x, cy - y, color);
    SetPixel(hdc, cx + x, cy - y, color);
    SetPixel(hdc, cx + y, cy + x, color);
    SetPixel(hdc, cx - y, cy + x, color);
    SetPixel(hdc, cx - y, cy - x, color);
    SetPixel(hdc, cx + y, cy - x, color);
}

// Midpoint Circle Algorithm (Bresenham variant)
void midpointCircle(HDC hdc, int cx, int cy, int radius, COLORREF color)
{
    int x = 0;
    int y = radius;
    int d = 1 - radius;

    draw8Points(hdc, cx, cy, x, y, color);

    while (x < y) {
        if (d < 0) {
            d += 2 * x + 3;
        } else {
            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
        draw8Points(hdc, cx, cy, x, y, color);
    }
}

// Direct Circle Method (using sqrt)
void circleDirectMethod(HDC hdc, int cx, int cy, int R, COLORREF color)
{
    double x = 0.0;
    double y = static_cast<double>(R);
    double R2 = static_cast<double>(R) * R;

    while (x < y) {
        draw8Points(hdc, cx, cy, static_cast<int>(round(x)), static_cast<int>(round(y)), color);
        x += 0.1;
        y = sqrt(R2 - x * x);
    }
}

// Polar Circle Method
void circlePolar(HDC hdc, int cx, int cy, int R, COLORREF color)
{
    int x = R;
    int y = 0;
    double theta = 0.0;
    double dtheta = 1.0 / R;

    draw8Points(hdc, cx, cy, x, y, color);

    while (x > y) {
        theta += dtheta;
        x = static_cast<int>(round(R * cos(theta)));
        y = static_cast<int>(round(R * sin(theta)));
        draw8Points(hdc, cx, cy, x, y, color);
    }
}

// Iterative Polar Circle Method
void circleIterative(HDC hdc, int cx, int cy, int R, COLORREF color)
{
    double x = static_cast<double>(R);
    double y = 0.0;
    double dtheta = 1.0 / R;
    double c = cos(dtheta);
    double s = sin(dtheta);

    draw8Points(hdc, cx, cy, static_cast<int>(round(x)), static_cast<int>(round(y)), color);

    while (x + 20.0 > y) {
        double xNew = (x * c) - (y * s);
        y = (x * s) + (y * c);
        x = xNew;
        draw8Points(hdc, cx, cy, static_cast<int>(round(x)), static_cast<int>(round(y)), color);
    }
}

// Modified Bresenham Circle (faster variant)
void circleFastBresenham(HDC hdc, int cx, int cy, int R, COLORREF color)
{
    int x = 0;
    int y = R;
    int d = 1 - R;
    int c1 = 3;
    int c2 = 5 - 2 * R;

    draw8Points(hdc, cx, cy, x, y, color);

    while (x < y) {
        if (d < 0) {
            d += c1;
            c2 += 2;
        } else {
            d += c2;
            c2 += 4;
            y--;
        }
        c1 += 2;
        x++;
        draw8Points(hdc, cx, cy, x, y, color);
    }
}

// Draws two symmetric points in a single quadrant
// quadrant = 0: top-right
// quadrant = 1: top-left
// quadrant = 2: bottom-left
// quadrant = 3: bottom-right
void drawQuadrant(HDC hdc, int cx, int cy, int x, int y, int quadrant, COLORREF color)
{
    switch (quadrant) {
    case 0:
        SetPixel(hdc, cx + x, cy - y, color);
        SetPixel(hdc, cx + y, cy - x, color);
        break;
    case 1:
        SetPixel(hdc, cx - x, cy - y, color);
        SetPixel(hdc, cx - y, cy - x, color);
        break;
    case 2:
        SetPixel(hdc, cx - x, cy + y, color);
        SetPixel(hdc, cx - y, cy + x, color);
        break;
    case 3:
        SetPixel(hdc, cx + x, cy + y, color);
        SetPixel(hdc, cx + y, cy + x, color);
        break;
    default:
        break;
    }
}

// Faster Bresenham circle restricted to a single quadrant
void circleFastBresenhamQuadrant(HDC hdc, int cx, int cy, int R, int quadrant, COLORREF color)
{
    int x = 0;
    int y = R;
    int d = 1 - R;
    int c1 = 3;
    int c2 = 5 - 2 * R;

    drawQuadrant(hdc, cx, cy, x, y, quadrant, color);

    while (x < y) {
        if (d < 0) {
            d += c1;
            c2 += 2;
        } else {
            d += c2;
            c2 += 4;
            y--;
        }
        c1 += 2;
        x++;
        drawQuadrant(hdc, cx, cy, x, y, quadrant, color);
    }
}

//-------------------------------------------------------------------------------------------------
// Ellipse-Drawing Algorithms
//-------------------------------------------------------------------------------------------------

// Draws four symmetric points around the ellipse center
void draw4Points(HDC hdc, int cx, int cy, int x, int y, COLORREF color)
{
    SetPixel(hdc, cx + x, cy + y, color);
    SetPixel(hdc, cx - x, cy + y, color);
    SetPixel(hdc, cx - x, cy - y, color);
    SetPixel(hdc, cx + x, cy - y, color);
}

// Direct Ellipse Method (using sqrt)
void ellipseDirect(HDC hdc, int cx, int cy, int a, int b, COLORREF color)
{
    int x = 0;
    int y = b;
    draw4Points(hdc, cx, cy, x, y, color);

    // Region 1: slope > -1
    while (x * (b * b) < (a * a) * y) {
        x++;
        y = static_cast<int>(round(sqrt((b * b * (a * a - x * x)) / (a * a))));
        draw4Points(hdc, cx, cy, x, y, color);
    }

    // Region 2: slope <= -1
    y = 0;
    x = a;
    draw4Points(hdc, cx, cy, x, y, color);

    while (x * (b * b) > (a * a) * y) {
        y++;
        x = static_cast<int>(round(sqrt((a * a * (b * b - y * y)) / (b * b))));
        draw4Points(hdc, cx, cy, x, y, color);
    }
}

// Polar Ellipse Method
void ellipsePolar(HDC hdc, int cx, int cy, int R1, int R2, COLORREF color)
{
    int x = R1;
    int y = 0;
    double theta = 0.0;
    double dtheta = 1.0 / R1;

    draw4Points(hdc, cx, cy, x, y, color);

    while (x + R1 > y) {
        theta += dtheta;
        x = static_cast<int>(round(R1 * cos(theta)));
        y = static_cast<int>(round(R2 * sin(theta)));
        draw4Points(hdc, cx, cy, x, y, color);
    }
}

//-------------------------------------------------------------------------------------------------
// Filling Algorithms
//-------------------------------------------------------------------------------------------------

// Fills a quarter-circle by drawing radial lines from center
// q = 1 → top-right, 2 → bottom-right, 3 → bottom-left, 4 → top-left
void fillWithLines(HDC hdc, int cx, int cy, int unusedA, int unusedB, int R, COLORREF color, int q)
{
    int x = R;
    int y = 0;
    double theta = 0.0;
    double dtheta = 1.0 / R;

    // Draw initial line segment in the chosen quadrant
    switch (q) {
    case 1: // top-right
        lineBresenham(hdc, cx, cy, cx + x, cy - y, color);
        break;
    case 2: // bottom-right
        lineBresenham(hdc, cx, cy, cx + x, cy + y, color);
        break;
    case 3: // bottom-left
        lineBresenham(hdc, cx, cy, cx - x, cy + y, color);
        break;
    case 4: // top-left
        lineBresenham(hdc, cx, cy, cx - x, cy - y, color);
        break;
    default:
        return; // invalid quadrant
    }

    while (x * R > y) {
        theta += dtheta;
        x = static_cast<int>(round(R * cos(theta)));
        y = static_cast<int>(round(R * sin(theta)));

        switch (q) {
        case 1: // top-right
            lineBresenham(hdc, cx, cy, cx + x, cy - y, color);
            break;
        case 2: // bottom-right
            lineBresenham(hdc, cx, cy, cx + x, cy + y, color);
            break;
        case 3: // bottom-left
            lineBresenham(hdc, cx, cy, cx - x, cy + y, color);
            break;
        case 4: // top-left
            lineBresenham(hdc, cx, cy, cx - x, cy - y, color);
            break;
        }
    }
}

// Fills a quarter-circle by drawing concentric arc segments
// q = 1 → top-right, 2 → bottom-right, 3 → bottom-left, 4 → top-left
void fillWithCircles(HDC hdc, int cx, int cy, int unusedA, int unusedB, int R, COLORREF color, int q)
{
    // Map q (1..4) to circleFastBresenhamQuadrant's indexing (0..3):
    //   q=1 (top-right)    → circleIndex = 0
    //   q=2 (bottom-right) → circleIndex = 3
    //   q=3 (bottom-left)  → circleIndex = 2
    //   q=4 (top-left)     → circleIndex = 1
    int circleIndex;
    switch (q) {
    case 1: circleIndex = 0; break;
    case 2: circleIndex = 3; break;
    case 3: circleIndex = 2; break;
    case 4: circleIndex = 1; break;
    default:
        return; // invalid quadrant
    }

    // Draw concentric “Fast Bresenham” arcs for each radius from 0 up to R-1
    for (int r = 0; r < R; r++) {
        circleFastBresenhamQuadrant(hdc, cx, cy, r, circleIndex, color);
    }
}

//-------------------------------------------------------------------------------------------------
// Point Clipping Inside a Circle
//-------------------------------------------------------------------------------------------------

// Colors a point blue if it lies inside the circle, otherwise red
void pointClip(HDC hdc, int x, int y, int radius, int centerX, int centerY, COLORREF color)
{
    double dx2 = pow(centerX - x, 2);
    double dy2 = pow(centerY - y, 2);
    double dist = sqrt(dx2 + dy2);

    if (dist <= radius) {
        SetPixel(hdc, x, y, RGB(0, 0, 255));
    } else {
        SetPixel(hdc, x, y, RGB(255, 0, 0));
    }
}

//-------------------------------------------------------------------------------------------------
// Window Procedure and Message Handling
//-------------------------------------------------------------------------------------------------

LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    static int currentShape = 0;                // ID of selected shape or action
    static COLORREF currentColor = RGB(0, 0, 0); // Current drawing color

    HDC hdc = GetDC(hwnd);

    switch (message) {
    case WM_COMMAND:
        // File → Save
        if (wParam == 1) {
            PerformSave("shapes_data.txt");
            ReleaseDC(hwnd, hdc);
            cout << "Action: Save shapes" << endl;
        }
        // File → Load
        else if (wParam == 2) {
            PerformLoad("shapes_data.txt");
            InvalidateRect(hwnd, NULL, TRUE);  // Trigger WM_PAINT
            cout << "Action: Load shapes" << endl;
        }
        // File → Clear
        else if (wParam == 3) {
            shapes.clear();
            InvalidateRect(hwnd, NULL, TRUE);
            cout << "Action: Clear canvas" << endl;
        }
        // Color Selection
        else if (wParam == 41) {
            currentColor = RGB(255, 0, 0);
            cout << "Action: Set color Red" << endl;
        }
        else if (wParam == 42) {
            currentColor = RGB(0, 0, 0);
            cout << "Action: Set color Black" << endl;
        }
        else if (wParam == 43) {
            currentColor = RGB(0, 0, 255);
            cout << "Action: Set color Blue" << endl;
        }
        else if (wParam == 44) {
            currentColor = RGB(255, 255, 255);
            cout << "Action: Set color White" << endl;
        }
        else if (wParam == 45) {
            currentColor = RGB(255, 63, 127);
            cout << "Action: Set color Custom Pink" << endl;
        }

        // Line Filling (IDs 101–104 map to quadrants 1–4)
        else if (wParam >= 101 && wParam <= 104) {
            int quadrant = (int)wParam - 100; // yields 1..4
            fillWithLines(hdc, xc, yc, xe, ye, R, currentColor, quadrant);
            cout << "Action: Fill quarter-circle with lines, Quadrant " << quadrant << endl;
        }
        // Circle Filling (IDs 201–204 map to quadrants 1–4)
        else if (wParam >= 201 && wParam <= 204) {
            int quadrant = (int)wParam - 200; // yields 1..4
            fillWithCircles(hdc, xc, yc, xe, ye, R, currentColor, quadrant);
            cout << "Action: Fill quarter-circle with circles, Quadrant " << quadrant << endl;
        }
        // Any other wParam is treated as selecting a shape/clipping mode
        else {
            currentShape = (int)wParam;
            switch (currentShape) {
            case 11: cout << "Mode: Draw DDA Line" << endl; break;
            case 12: cout << "Mode: Draw Midpoint Line" << endl; break;
            case 13: cout << "Mode: Draw Parametric Line" << endl; break;
            case 14: cout << "Mode: Draw Circle (Direct)" << endl; break;
            case 15: cout << "Mode: Draw Circle (Polar)" << endl; break;
            case 16: cout << "Mode: Draw Circle (Iterative)" << endl; break;
            case 17: cout << "Mode: Draw Circle (Midpoint)" << endl; break;
            case 18: cout << "Mode: Draw Circle (Fast Bresenham)" << endl; break;
            case 19: cout << "Mode: Draw Ellipse (Direct)" << endl; break;
            case 20: cout << "Mode: Draw Ellipse (Polar)" << endl; break;
            case 51: cout << "Mode: Define Clipping Window" << endl; break;
            case 52: cout << "Mode: Clip Line in Window" << endl; break;
            case 53: cout << "Mode: Clip Point in Window" << endl; break;
            default: break;
            }
        }
        break;

    case WM_LBUTTONDOWN:
        xc = LOWORD(lParam);
        yc = HIWORD(lParam);

        // If Clipping Window mode is selected (ID = 51)
        if (currentShape == 51) {
            windowX = xc;
            windowY = yc;
            cout << "Action: Set clipping window center (" << xc << ", " << yc << ")" << endl;
        }
        // If Point Clipping mode is selected (ID = 53)
        else if (currentShape == 53) {
            pointClip(hdc, xc, yc, windowR, windowX, windowY, currentColor);
            cout << "Action: Clip point at (" << xc << ", " << yc << ")" << endl;
        }
        break;

    case WM_RBUTTONDOWN:
        xe = LOWORD(lParam);
        ye = HIWORD(lParam);

        // Drawing operations based on selected shape ID

        // DDALine (ID = 11)
        if (currentShape == 11) {
            DDALine(hdc, xc, yc, xe, ye, currentColor);
            shapes.push_back({ "DDALine", xc, yc, xe, ye, 0, 0, currentColor });
            cout << "Drawn: DDA Line from (" << xc << "," << yc << ") to (" << xe << "," << ye << ")" << endl;
        }
        // Midpoint Line (ID = 12)
        else if (currentShape == 12) {
            lineBresenham(hdc, xc, yc, xe, ye, currentColor);
            shapes.push_back({ "MidLine", xc, yc, xe, ye, 0, 0, currentColor });
            cout << "Drawn: Midpoint Line from (" << xc << "," << yc << ") to (" << xe << "," << ye << ")" << endl;
        }
        // Parametric Line (ID = 13)
        else if (currentShape == 13) {
            parametricLine(hdc, xc, yc, xe, ye, 0, 0, 0, currentColor);
            shapes.push_back({ "ParamLine", xc, yc, xe, ye, 0, 0, currentColor });
            cout << "Drawn: Parametric Line from (" << xc << "," << yc << ") to (" << xe << "," << ye << ")" << endl;
        }

        // Circle Direct (ID = 14)
        else if (currentShape == 14) {
            R = static_cast<int>(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circleDirectMethod(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleDirect", xc, yc, 0, 0, R, 0, currentColor });
            cout << "Drawn: Circle (Direct) center (" << xc << "," << yc << "), radius " << R << endl;
        }
        // Circle Polar (ID = 15)
        else if (currentShape == 15) {
            R = static_cast<int>(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circlePolar(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CirclePolar", xc, yc, 0, 0, R, 0, currentColor });
            cout << "Drawn: Circle (Polar) center (" << xc << "," << yc << "), radius " << R << endl;
        }
        // Circle Iterative (ID = 16)
        else if (currentShape == 16) {
            R = static_cast<int>(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circleIterative(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleIter", xc, yc, 0, 0, R, 0, currentColor });
            cout << "Drawn: Circle (Iterative) center (" << xc << "," << yc << "), radius " << R << endl;
        }
        // Circle Midpoint (ID = 17)
        else if (currentShape == 17) {
            R = static_cast<int>(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            midpointCircle(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleMid", xc, yc, 0, 0, R, 0, currentColor });
            cout << "Drawn: Circle (Midpoint) center (" << xc << "," << yc << "), radius " << R << endl;
        }
        // Circle Fast Bresenham (ID = 18)
        else if (currentShape == 18) {
            R = static_cast<int>(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circleFastBresenham(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleFast", xc, yc, 0, 0, R, 0, currentColor });
            cout << "Drawn: Circle (Fast Bresenham) center (" << xc << "," << yc << "), radius " << R << endl;
        }

        // Ellipse Direct (ID = 19)
        else if (currentShape == 19) {
            int a = abs(xe - xc);
            int b = abs(ye - yc);
            ellipseDirect(hdc, xc, yc, a, b, currentColor);
            shapes.push_back({ "EllipseDirect", xc, yc, 0, 0, a, b, currentColor });
            cout << "Drawn: Ellipse (Direct) center (" << xc << "," << yc << "), axes = (" << a << ", " << b << ")" << endl;
        }

        // Clipping Window Definition (ID = 51)
        else if (currentShape == 51) {
            windowX = xc;
            windowY = yc;
            R = static_cast<int>(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            windowR = R;
            circlePolar(hdc, xc, yc, R, currentColor);
            cout << "Defined clipping window center (" << xc << "," << yc << "), radius " << R << endl;
        }
        // Line Clipping (within window) (ID = 52)
        else if (currentShape == 52) {
            parametricLine(hdc, xc, yc, xe, ye, windowX, windowY, windowR, currentColor);
            cout << "Clipped line from (" << xc << "," << yc << ") to (" << xe << "," << ye << ") within window" << endl;
        }

        break;

    case WM_RBUTTONUP:
        // Polar Ellipse requires two radii (ID = 20)
        if (currentShape == 20) {
            xe2 = LOWORD(lParam);
            ye2 = HIWORD(lParam);
            R  = static_cast<int>(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            R2 = static_cast<int>(sqrt(pow(xe2 - xc, 2) + pow(ye2 - yc, 2)));
            ellipsePolar(hdc, xc, yc, R, R2, currentColor);
            shapes.push_back({ "EllipsePolar", xc, yc, 0, 0, R, R2, currentColor });
            cout << "Drawn: Ellipse (Polar) center (" << xc << "," << yc << "), radii = (" << R << ", " << R2 << ")" << endl;
        }
        break;

    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdcPaint = BeginPaint(hwnd, &ps);

        // Redraw all stored shapes
        for (const auto& s : shapes) {
            if (s.type == "DDALine") {
                DDALine(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "MidLine") {
                lineBresenham(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "ParamLine") {
                parametricLine(hdcPaint, s.x1, s.y1, s.x2, s.y2, 0, 0, 0, s.color);
            }
            else if (s.type == "CircleDirect") {
                circleDirectMethod(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CirclePolar") {
                circlePolar(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CircleIter") {
                circleIterative(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CircleMid") {
                midpointCircle(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CircleFast") {
                circleFastBresenham(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "EllipseDirect") {
                ellipseDirect(hdcPaint, s.x1, s.y1, s.a, s.b, s.color);
            }
            else if (s.type == "EllipsePolar") {
                ellipsePolar(hdcPaint, s.x1, s.y1, s.a, s.b, s.color);
            }
        }

        EndPaint(hwnd, &ps);
        break;
    }

    case WM_CREATE:
        // Create and attach the menu to the window
        menu(hwnd);
        break;

    case WM_DESTROY:
        PostQuitMessage(0);
        break;

    default:
        return DefWindowProc(hwnd, message, wParam, lParam);
    }

    return 0;
}

//-------------------------------------------------------------------------------------------------
// WinMain: Application Entry Point
//-------------------------------------------------------------------------------------------------

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    // Allocate a console so that we can print messages to it
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    ShowWindow(GetConsoleWindow(), SW_SHOW);

    // Register the window class
    WNDCLASSEX wincl = {};
    wincl.cbSize        = sizeof(WNDCLASSEX);
    wincl.style         = CS_DBLCLKS;
    wincl.lpfnWndProc   = WindowProcedure;
    wincl.hInstance     = hInstance;
    wincl.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wincl.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);
    wincl.hCursor       = LoadCursor(NULL, IDC_CROSS);
    wincl.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
    wincl.lpszClassName = _T("CodeBlocksWindowsApp");
    wincl.lpszMenuName  = NULL;

    if (!RegisterClassEx(&wincl)) {
        return 0;
    }

    // Create the application window
    HWND hwnd = CreateWindowEx(
        0,
        _T("CodeBlocksWindowsApp"),
        _T("Computer Graphics Project 2025"),
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT,
        CW_USEDEFAULT,
        544,
        375,
        HWND_DESKTOP,
        NULL,
        hInstance,
        NULL
    );

    ShowWindow(hwnd, nCmdShow);

    // Message loop
    MSG messages;
    while (GetMessage(&messages, NULL, 0, 0)) {
        TranslateMessage(&messages);
        DispatchMessage(&messages);
    }

    return messages.wParam;
}

//-------------------------------------------------------------------------------------------------
// Menu Creation
//-------------------------------------------------------------------------------------------------

void menu(HWND hwnd)
{
    hMenu = CreateMenu();

    // File Menu: Save, Load, Clear
    HMENU hFileMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFileMenu, "File");
    AppendMenu(hFileMenu, MF_STRING, 1, "Save");
    AppendMenu(hFileMenu, MF_STRING, 2, "Load");
    AppendMenu(hFileMenu, MF_STRING, 3, "Clear");

    // Line Menu: DDA, Midpoint, Parametric
    HMENU hLineMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hLineMenu, "Line");
    AppendMenu(hLineMenu, MF_STRING, 11, "DDA");
    AppendMenu(hLineMenu, MF_STRING, 12, "Midpoint");
    AppendMenu(hLineMenu, MF_STRING, 13, "Parametric");

    // Circle Menu: Direct, Polar, Iterative Polar, Midpoint, Modified Midpoint
    HMENU hCircleMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCircleMenu, "Circle");
    AppendMenu(hCircleMenu, MF_STRING, 14, "Direct");
    AppendMenu(hCircleMenu, MF_STRING, 15, "Polar");
    AppendMenu(hCircleMenu, MF_STRING, 16, "Iterative Polar");
    AppendMenu(hCircleMenu, MF_STRING, 17, "Midpoint");
    AppendMenu(hCircleMenu, MF_STRING, 18, "Modified Midpoint");

    // Ellipse Menu: Direct, Polar
    HMENU hEllipseMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hEllipseMenu, "Ellipse");
    AppendMenu(hEllipseMenu, MF_STRING, 19, "Direct");
    AppendMenu(hEllipseMenu, MF_STRING, 20, "Polar");

    // Color Menu: choose drawing color
    HMENU hColorMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hColorMenu, "Color");
    AppendMenu(hColorMenu, MF_STRING, 41, "Red");
    AppendMenu(hColorMenu, MF_STRING, 42, "Black");
    AppendMenu(hColorMenu, MF_STRING, 43, "Blue");
    AppendMenu(hColorMenu, MF_STRING, 44, "White");
    AppendMenu(hColorMenu, MF_STRING, 45, "Custom Pink");

    // Clipping Menu: Window, Line, Point
    HMENU hClippingMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClippingMenu, "Clipping");
    AppendMenu(hClippingMenu, MF_STRING, 51, "Window");
    AppendMenu(hClippingMenu, MF_STRING, 52, "Line");
    AppendMenu(hClippingMenu, MF_STRING, 53, "Point");

    // Line Filling Menu (IDs 101–104 correspond to quadrants 1–4)
    HMENU hLineFillingMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hLineFillingMenu, "Line Filling");
    AppendMenu(hLineFillingMenu, MF_STRING, 101, "Quadrant 1");
    AppendMenu(hLineFillingMenu, MF_STRING, 102, "Quadrant 2");
    AppendMenu(hLineFillingMenu, MF_STRING, 103, "Quadrant 3");
    AppendMenu(hLineFillingMenu, MF_STRING, 104, "Quadrant 4");

    // Circle Filling Menu (IDs 201–204 correspond to quadrants 1–4)
    HMENU hCircleFillingMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCircleFillingMenu, "Circle Filling");
    AppendMenu(hCircleFillingMenu, MF_STRING, 201, "Quadrant 1");
    AppendMenu(hCircleFillingMenu, MF_STRING, 202, "Quadrant 2");
    AppendMenu(hCircleFillingMenu, MF_STRING, 203, "Quadrant 3");
    AppendMenu(hCircleFillingMenu, MF_STRING, 204, "Quadrant 4");

    SetMenu(hwnd, hMenu);
}

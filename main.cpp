#include <tchar.h>
#include <windows.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stack>
#include <algorithm>

using namespace std;

//-----------------------------------------------------------------------------
// Requirement 1 & 2 & 3 & 4: Window & Console Setup
// 1. White background is set in WNDCLASSEX.hbrBackground.
// 2. Crosshair cursor used.
// 3. Mouse-only interaction (no keyboard handlers).
// 4. AllocConsole + freopen for console I/O.
//-----------------------------------------------------------------------------

LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);

//-----------------------------------------------------------------------------
// Data Structures and Global Variables
//-----------------------------------------------------------------------------

struct Shape {
    string type;          // e.g., "DDA", "CirclePolar", "ClipCircleLine", etc.
    int x1, y1;           // primary coordinates
    int x2, y2;           // secondary coordinates
    int a, b;             // extra parameters (radius, etc.)
    int c;                // auxiliary
    COLORREF color;       // drawing color
    vector<POINT> verts;  // for polygons or splines
};

static vector<Shape> shapes;

// Mouse click coordinates
static int xc, yc;      // first click
static int xe, ye;      // second click
static int xe2, ye2;    // third click (for polar ellipse)
static int R, R2;       // radii or aux

// Clipping windows
static int windowX, windowY, windowR;              // circular clipping
static int rectXmin, rectYmin, rectXmax, rectYmax; // rectangular clipping
static int sqXmin, sqYmin, sqSide;                  // square clipping

static vector<POINT> currentVerts;        // for polygon or spline
static HMENU hMenu;

//-----------------------------------------------------------------------------
// Mode IDs (named constants instead of raw numbers)
//-----------------------------------------------------------------------------

// File / Clear
const int ID_SAVE                = 1;
const int ID_LOAD                = 2;
const int ID_CLEAR               = 3;

// Color submenu
const int ID_COLOR_RED           = 41;
const int ID_COLOR_BLACK         = 42;
const int ID_COLOR_BLUE          = 43;
const int ID_COLOR_WHITE         = 44;
const int ID_COLOR_PINK          = 45;

// Flood fill
const int ID_FLOOD_REC           = 70;
const int ID_FLOOD_STACK         = 71;

// Bezier rectangle fill
const int ID_FILL_BEZIERRECT     = 80;

// Polygon fill
const int ID_FILL_CONVEX_POLY    = 90;
const int ID_FILL_NONCONVEX_POLY = 91;

// Cardinal spline
const int ID_CARDINAL_SPLINE     = 95;

// Rectangular clipping
const int ID_CLIP_RECT_WINDOW    = 100;
const int ID_CLIP_RECT_POINT     = 101;
const int ID_CLIP_RECT_LINE      = 102;
const int ID_CLIP_RECT_POLY      = 103;

// Square clipping
const int ID_CLIP_SQUARE_WINDOW  = 110;
const int ID_CLIP_SQUARE_POINT   = 111;
const int ID_CLIP_SQUARE_LINE    = 112;

// Line Filling (unused in this request)
const int ID_FILL_LINE_Q1        = 210;
const int ID_FILL_LINE_Q2        = 211;
const int ID_FILL_LINE_Q3        = 212;
const int ID_FILL_LINE_Q4        = 213;

// Circle Filling (lines)
const int ID_FILL_CIRCLE_Q1      = 201;
const int ID_FILL_CIRCLE_Q2      = 202;
const int ID_FILL_CIRCLE_Q3      = 203;
const int ID_FILL_CIRCLE_Q4      = 204;

// Line algorithms
const int ID_LINE_DDA            = 11;
const int ID_LINE_MID            = 12;
const int ID_LINE_PARAM          = 13;

// Circle algorithms
const int ID_CIRCLE_DIRECT       = 14;
const int ID_CIRCLE_POLAR        = 15;
const int ID_CIRCLE_ITERATIVE    = 16;
const int ID_CIRCLE_MIDPOINT     = 17;
const int ID_CIRCLE_FAST         = 18;

// Hermite square fill
const int ID_HERMITE_SQUARE      = 60;

// Ellipse algorithms
const int ID_ELLIPSE_DIRECT      = 19;
const int ID_ELLIPSE_POLAR       = 20;
const int ID_ELLIPSE_MIDPOINT    = 21;

// Circular clipping (bonus)
const int ID_CLIP_CIRCLE_WINDOW  = 51;
const int ID_CLIP_CIRCLE_LINE    = 52;
const int ID_CLIP_CIRCLE_POINT   = 53;

//-----------------------------------------------------------------------------
// Utility Functions: LERP and Bezier
//-----------------------------------------------------------------------------

static double LERP(double s, double e, double t) {
    return (1.0 - t) * s + t * e;
}

static double quadraticBezier(double x, double y, double z, double t) {
    double a = LERP(x, y, t);
    double b = LERP(y, z, t);
    return LERP(a, b, t);
}

//-----------------------------------------------------------------------------
// Requirement 9: Line Algorithms [DDA, Midpoint, Parametric]
//-----------------------------------------------------------------------------

static void lineBresenham(HDC hdc, int x0, int y0, int x1, int y1, COLORREF color) {
    int deltaX = abs(x1 - x0);
    int deltaY = abs(y1 - y0);
    int stepX = (x0 < x1) ? 1 : -1;
    int stepY = (y0 < y1) ? 1 : -1;
    int decision = deltaX - deltaY;

    while (true) {
        SetPixel(hdc, x0, y0, color);  // Plot current point
        if (x0 == x1 && y0 == y1) {
            break;  // Reached the end point
        }
        int doubledDecision = 2 * decision;
        if (doubledDecision > -deltaY) {
            decision -= deltaY;
            x0 += stepX;
        }
        if (doubledDecision < deltaX) {
            decision += deltaX;
            y0 += stepY;
        }
    }
}
static void DDALine(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    int dx = x2 - x1, dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    if (steps == 0) {
        SetPixel(hdc, x1, y1, color);
        return;
    }
    float Xinc = dx / (float)steps;
    float Yinc = dy / (float)steps;
    float X = (float)x1, Y = (float)y1;
    for (int i = 0; i <= steps; i++) {
        SetPixel(hdc, (int)round(X), (int)round(Y), color);
        X += Xinc;
        Y += Yinc;
    }
}

static void parametricLine(HDC hdc, int x1, int y1, int x2, int y2, int cx, int cy, int radius, COLORREF color) {
    int dx = x2 - x1, dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    if (steps == 0) {
        if (!radius || ((x1 - cx)*(x1 - cx) + (y1 - cy)*(y1 - cy) <= radius*radius)) {
            SetPixel(hdc, x1, y1, color);
        }
        return;
    }
    double tStep = 1.0 / steps;
    for (double t = 0.0; t <= 1.0; t += tStep) {
        int x = (int)round(LERP(x1, x2, t));
        int y = (int)round(LERP(y1, y2, t));
        if (!radius || ((x - cx)*(x - cx) + (y - cy)*(y - cy) <= radius*radius)) {
            SetPixel(hdc, x, y, color);
        }
    }
}

//-----------------------------------------------------------------------------
// Requirement 10: Circle Algorithms [Direct, Polar, Iterative Polar, Midpoint, Fast Bresenham]
//-----------------------------------------------------------------------------

static void draw8Points(HDC hdc, int cx, int cy, int x, int y, COLORREF color) {
    SetPixel(hdc, cx + x, cy + y, color);
    SetPixel(hdc, cx - x, cy + y, color);
    SetPixel(hdc, cx - x, cy - y, color);
    SetPixel(hdc, cx + x, cy - y, color);
    SetPixel(hdc, cx + y, cy + x, color);
    SetPixel(hdc, cx - y, cy + x, color);
    SetPixel(hdc, cx - y, cy - x, color);
    SetPixel(hdc, cx + y, cy - x, color);
}

static void midpointCircle(HDC hdc, int cx, int cy, int radius, COLORREF color) {
    int x = 0, y = radius;
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

static void circleDirectMethod(HDC hdc, int cx, int cy, int R, COLORREF color) {
    double x = 0.0, y = (double)R;
    double R2 = (double)R * R;
    while (x < y) {
        draw8Points(hdc, cx, cy, (int)round(x), (int)round(y), color);
        x += 0.1;
        y = sqrt(R2 - x * x);
        draw8Points(hdc, cx, cy, (int)round(x), (int)round(y), color);
    }
}

static void circlePolar(HDC hdc, int cx, int cy, int R, COLORREF color) {
    int x = R, y = 0;
    double theta = 0.0, dtheta = 1.0 / R;
    draw8Points(hdc, cx, cy, x, y, color);
    while (x > y) {
        theta += dtheta;
        x = (int)round(R * cos(theta));
        y = (int)round(R * sin(theta));
        draw8Points(hdc, cx, cy, x, y, color);
    }
}

static void circleIterative(HDC hdc, int cx, int cy, int R, COLORREF color) {
    double x = (double)R, y = 0.0;
    double dtheta = 1.0 / R;
    double c = cos(dtheta), s = sin(dtheta);
    draw8Points(hdc, cx, cy, (int)round(x), (int)round(y), color);
    while (x + 20.0 > y) {
        double xNew = x * c - y * s;
        y = x * s + y * c;
        x = xNew;
        draw8Points(hdc, cx, cy, (int)round(x), (int)round(y), color);
    }
}

static void circleFastBresenham(HDC hdc, int cx, int cy, int R, COLORREF color) {
    int x = 0, y = R;
    int d = 1 - R, c1 = 3, c2 = 5 - 2 * R;
    draw8Points(hdc, cx, cy, x, y, color);
    while (x < y) {
        if (d < 0) {
            d += c1; c2 += 2;
        } else {
            d += c2; c2 += 4; y--;
        }
        x++; c1 += 2;
        draw8Points(hdc, cx, cy, x, y, color);
    }
}

//-----------------------------------------------------------------------------
// Requirement 11: Filling Circle with Lines (quarter selection)
//-----------------------------------------------------------------------------

static void drawQuadrant(HDC hdc, int cx, int cy, int x, int y, int quadrant, COLORREF color) {
    switch (quadrant) {
        case 1:
            SetPixel(hdc, cx + x, cy - y, color);
            SetPixel(hdc, cx + y, cy - x, color);
            break;
        case 2:
            SetPixel(hdc, cx - x, cy - y, color);
            SetPixel(hdc, cx - y, cy - x, color);
            break;
        case 3:
            SetPixel(hdc, cx - x, cy + y, color);
            SetPixel(hdc, cx - y, cy + x, color);
            break;
        case 4:
            SetPixel(hdc, cx + x, cy + y, color);
            SetPixel(hdc, cx + y, cy + x, color);
            break;
    }
}

static void fillWithLines(HDC hdc, int cx, int cy, int unusedA, int unusedB, int R, COLORREF color, int quad) {
    int x = R, y = 0;
    double theta = 0.0, dtheta = 1.0 / R;
    switch (quad) {
        case 1: lineBresenham(hdc, cx, cy, cx + x, cy - y, color); break;
        case 2: lineBresenham(hdc, cx, cy, cx + x, cy + y, color); break;
        case 3: lineBresenham(hdc, cx, cy, cx - x, cy + y, color); break;
        case 4: lineBresenham(hdc, cx, cy, cx - x, cy - y, color); break;
    }
    while (x * R > y) {
        theta += dtheta;
        x = (int)round(R * cos(theta));
        y = (int)round(R * sin(theta));
        switch (quad) {
            case 1: lineBresenham(hdc, cx, cy, cx + x, cy - y, color); break;
            case 2: lineBresenham(hdc, cx, cy, cx + x, cy + y, color); break;
            case 3: lineBresenham(hdc, cx, cy, cx - x, cy + y, color); break;
            case 4: lineBresenham(hdc, cx, cy, cx - x, cy - y, color); break;
        }
    }
}

//-----------------------------------------------------------------------------
// Requirement 12: Filling Circle with Smaller Circles (quarter selection)
//-----------------------------------------------------------------------------

static void fillWithCircles(HDC hdc, int cx, int cy, int unusedA, int unusedB, int R, COLORREF color, int quad) {
    for (int r = 0; r < R; r++) {
        drawQuadrant(hdc, cx, cy, 0, r, quad, color);
        int x = 0, y = r;
        int d = 1 - r, c1 = 3, c2 = 5 - 2 * r;
        while (x < y) {
            if (d < 0) {
                d += c1; c2 += 2;
            } else {
                d += c2; c2 += 4; y--;
            }
            x++; c1 += 2;
            drawQuadrant(hdc, cx, cy, x, y, quad, color);
        }
    }
}

//-----------------------------------------------------------------------------
// Requirement 13: Filling Square with Hermite Curve (vertical orientation)
//-----------------------------------------------------------------------------

static POINT EvaluateHermite(float t, POINT p0, POINT p1, POINT r0, POINT r1) {
    float h1 =  2*t*t*t - 3*t*t + 1;
    float h2 = -2*t*t*t + 3*t*t;
    float h3 =  t*t*t - 2*t*t + t;
    float h4 =  t*t*t - t*t;
    POINT p;
    p.x = (int)round(h1*p0.x + h2*p1.x + h3*r0.x + h4*r1.x);
    p.y = (int)round(h1*p0.y + h2*p1.y + h3*r0.y + h4*r1.y);
    return p;
}

static void FillHermiteSquare(HDC hdc, POINT topLeft, int size, COLORREF color) {
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    const int COUNT = 50;
    for (int i = 0; i <= COUNT; i++) {
        float t = (float)i / COUNT;
        int x = topLeft.x + (int)round(t*size);
        POINT p0 = { x, topLeft.y };
        POINT p1 = { x, topLeft.y + size };
        POINT r0 = { (int)round(50.0f * sinf(3.1415f * t)), 100 };
        POINT r1 = { (int)round(-50.0f * sinf(3.1415f * t)), -100 };

        const int STEPS = 100;
        POINT prev = EvaluateHermite(0.0f, p0, p1, r0, r1);
        for (int j = 1; j <= STEPS; j++) {
            float tj = (float)j / STEPS;
            POINT curr = EvaluateHermite(tj, p0, p1, r0, r1);
            MoveToEx(hdc, prev.x, prev.y, NULL);
            LineTo(hdc, curr.x, curr.y);
            prev = curr;
        }
    }

    // Square outline
    MoveToEx(hdc, topLeft.x, topLeft.y, NULL);
    LineTo(hdc, topLeft.x + size, topLeft.y);
    LineTo(hdc, topLeft.x + size, topLeft.y + size);
    LineTo(hdc, topLeft.x, topLeft.y + size);
    LineTo(hdc, topLeft.x, topLeft.y);

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

//-----------------------------------------------------------------------------
// Requirement 14: Filling Rectangle with Horizontal Bezier Curve
//-----------------------------------------------------------------------------

static void FillBezierRectangle(HDC hdc, POINT topLeft, POINT bottomRight, COLORREF color) {
    int x1 = topLeft.x, y1 = topLeft.y;
    int x2 = bottomRight.x, y2 = bottomRight.y;
    if (x2 < x1) swap(x1, x2);
    if (y2 < y1) swap(y1, y2);

    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    int cx = (x1 + x2) / 2;
    int cy = (y1 + y2) / 2;
    const int SEGMENTS = 50;

    for (int y = y1; y <= y2; y++) {
        double px0 = x1, py0 = y;
        double px1 = cx, py1 = cy;
        double px2 = x2, py2 = y;
        double prevX = px0, prevY = py0;
        for (int i = 1; i <= SEGMENTS; i++) {
            double t = (double)i / SEGMENTS;
            double xt = quadraticBezier(px0, px1, px2, t);
            double yt = quadraticBezier(py0, py1, py2, t);
            MoveToEx(hdc, (int)round(prevX), (int)round(prevY), NULL);
            LineTo(hdc, (int)round(xt), (int)round(yt));
            prevX = xt;
            prevY = yt;
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

//-----------------------------------------------------------------------------
// Requirement 15: Convex & Non-Convex Polygon Filling
//-----------------------------------------------------------------------------

static void fillConvexPolygon(HDC hdc, const vector<POINT>& verts, COLORREF color) {
    if (verts.size() < 3) return;
    int ymin = (int)verts[0].y, ymax = (int)verts[0].y;
    for (auto &p : verts) {
        int py = (int)p.y;
        ymin = min(ymin, py);
        ymax = max(ymax, py);
    }
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    for (int y = ymin; y <= ymax; y++) {
        double leftX = 1e9, rightX = -1e9;
        for (size_t i = 0; i < verts.size(); i++) {
            POINT p1 = verts[i];
            POINT p2 = verts[(i + 1) % verts.size()];
            int p1y = (int)p1.y, p2y = (int)p2.y;
            if (p1y == p2y) continue;
            int edgeYmin = min(p1y, p2y);
            int edgeYmax = max(p1y, p2y);
            if (y < edgeYmin || y > edgeYmax) continue;
            double x = p1.x + (double)(y - p1.y) * (double)(p2.x - p1.x) / (double)(p2.y - p1.y);
            leftX  = min(leftX,  x);
            rightX = max(rightX, x);
        }
        if (leftX <= rightX) {
            MoveToEx(hdc, (int)ceil(leftX), y, NULL);
            LineTo(hdc, (int)floor(rightX), y);
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

static void fillNonConvexPolygon(HDC hdc, const vector<POINT>& verts, COLORREF color) {
    if (verts.size() < 3) return;
    int ymin = (int)verts[0].y, ymax = (int)verts[0].y;
    for (auto &p : verts) {
        int py = (int)p.y;
        ymin = min(ymin, py);
        ymax = max(ymax, py);
    }
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    for (int y = ymin; y <= ymax; y++) {
        vector<double> intersections;
        for (size_t i = 0; i < verts.size(); i++) {
            POINT p1 = verts[i];
            POINT p2 = verts[(i + 1) % verts.size()];
            int p1y = (int)p1.y, p2y = (int)p2.y;
            if (p1y == p2y) continue;
            int edgeYmin = min(p1y, p2y);
            int edgeYmax = max(p1y, p2y);
            if (y < edgeYmin || y >= edgeYmax) continue;
            double x = p1.x + (double)(y - p1.y) * (double)(p2.x - p1.x) / (double)(p2.y - p1.y);
            intersections.push_back(x);
        }
        sort(intersections.begin(), intersections.end());
        for (size_t i = 0; i + 1 < intersections.size(); i += 2) {
            int xStart = (int)ceil(intersections[i]);
            int xEnd   = (int)floor(intersections[i + 1]);
            if (xStart <= xEnd) {
                MoveToEx(hdc, xStart, y, NULL);
                LineTo(hdc, xEnd, y);
            }
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

//-----------------------------------------------------------------------------
// Requirement 16: Recursive & Non-Recursive Flood Fill
//-----------------------------------------------------------------------------

static void floodFillRec(HDC hdc, int x, int y, COLORREF targetColor, COLORREF fillColor) {
    COLORREF curr = GetPixel(hdc, x, y);
    if (curr != targetColor || curr == fillColor) return;
    SetPixel(hdc, x, y, fillColor);
    floodFillRec(hdc, x + 1, y, targetColor, fillColor);
    floodFillRec(hdc, x - 1, y, targetColor, fillColor);
    floodFillRec(hdc, x, y + 1, targetColor, fillColor);
    floodFillRec(hdc, x, y - 1, targetColor, fillColor);
}

static void floodFillStack(HDC hdc, int x, int y, COLORREF targetColor, COLORREF fillColor) {
    if (targetColor == fillColor) return;
    stack<POINT> stk;
    stk.push({ x, y });
    while (!stk.empty()) {
        POINT p = stk.top(); stk.pop();
        COLORREF curr = GetPixel(hdc, p.x, p.y);
        if (curr != targetColor || curr == fillColor) continue;
        SetPixel(hdc, p.x, p.y, fillColor);
        stk.push({ p.x + 1, p.y });
        stk.push({ p.x - 1, p.y });
        stk.push({ p.x, p.y + 1 });
        stk.push({ p.x, p.y - 1 });
    }
}

//-----------------------------------------------------------------------------
// Requirement 17: Cardinal Spline Curve
//-----------------------------------------------------------------------------

static void drawCardinalSpline(HDC hdc, const vector<POINT>& pts, COLORREF color) {
    size_t n = pts.size();
    if (n < 3) return;
    vector<POINT> tangents(n);
    tangents[0] = { 0, 0 };
    tangents[n - 1] = { 0, 0 };
    for (size_t i = 1; i < n - 1; i++) {
        tangents[i].x = (int)round(0.5 * (pts[i + 1].x - pts[i - 1].x));
        tangents[i].y = (int)round(0.5 * (pts[i + 1].y - pts[i - 1].y));
    }

    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    const int STEPS_SEG = 50;
    for (size_t i = 0; i < n - 1; i++) {
        POINT P0 = pts[i];
        POINT P1 = pts[i + 1];
        POINT R0 = tangents[i];
        POINT R1 = tangents[i + 1];

        POINT prev = EvaluateHermite(0.0f, P0, P1, R0, R1);
        for (int s = 1; s <= STEPS_SEG; s++) {
            float t = (float)s / STEPS_SEG;
            POINT curr = EvaluateHermite(t, P0, P1, R0, R1);
            MoveToEx(hdc, prev.x, prev.y, NULL);
            LineTo(hdc, curr.x, curr.y);
            prev = curr;
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

//-----------------------------------------------------------------------------
// Requirement 18: Ellipse Algorithms [Direct, Polar, Midpoint]
//-----------------------------------------------------------------------------

static void draw4Points(HDC hdc, int cx, int cy, int x, int y, COLORREF color) {
    SetPixel(hdc, cx + x, cy + y, color);
    SetPixel(hdc, cx - x, cy + y, color);
    SetPixel(hdc, cx - x, cy - y, color);
    SetPixel(hdc, cx + x, cy - y, color);
}

static void ellipseDirect(HDC hdc, int cx, int cy, int a, int b, COLORREF color) {
    int x = 0, y = b;
    draw4Points(hdc, cx, cy, x, y, color);
    while (x * (b * b) < (a * a) * y) {
        x++;
        y = (int)round(sqrt((double)(b * b * (a * a - x * x)) / (double)(a * a)));
        draw4Points(hdc, cx, cy, x, y, color);
    }
    y = 0;
    x = a;
    draw4Points(hdc, cx, cy, x, y, color);
    while (x * (b * b) > (a * a) * y) {
        y++;
        x = (int)round(sqrt((double)(a * a * (b * b - y * y)) / (double)(b * b)));
        draw4Points(hdc, cx, cy, x, y, color);
    }
}

static void ellipsePolar(HDC hdc, int cx, int cy, int R1, int R2, COLORREF color) {
    int x = R1, y = 0;
    double theta = 0.0, dtheta = 1.0 / R1;
    draw4Points(hdc, cx, cy, x, y, color);
    while (x + R1 > y) {
        theta += dtheta;
        x = (int)round(R1 * cos(theta));
        y = (int)round(R2 * sin(theta));
        draw4Points(hdc, cx, cy, x, y, color);
    }
}

static void ellipseMidpoint(HDC hdc, int cx, int cy, int a, int b, COLORREF color) {
    int x = 0, y = b;
    double a2 = (double)a * a, b2 = (double)b * b;
    double d1 = b2 - (a2 * b) + (0.25 * a2);
    draw4Points(hdc, cx, cy, x, y, color);
    double dx = 2 * b2 * x;
    double dy = 2 * a2 * y;

    // Region 1
    while (dx < dy) {
        x++;
        dx += 2 * b2;
        if (d1 < 0) {
            d1 += b2 + dx;
        } else {
            y--;
            dy -= 2 * a2;
            d1 += b2 + dx - dy;
        }
        draw4Points(hdc, cx, cy, x, y, color);
    }
    // Region 2
    double d2 = b2 * (x + 0.5) * (x + 0.5) + a2 * (y - 1) * (y - 1) - a2 * b2;
    while (y > 0) {
        y--;
        dy -= 2 * a2;
        if (d2 > 0) {
            d2 += a2 - dy;
        } else {
            x++;
            dx += 2 * b2;
            d2 += a2 - dy + dx;
        }
        draw4Points(hdc, cx, cy, x, y, color);
    }
}

//-----------------------------------------------------------------------------
// Requirement 19: Rectangular Clipping [Point, Line, Polygon]
//-----------------------------------------------------------------------------

static void drawRectangleWindow(HDC hdc, int xmin, int ymin, int xmax, int ymax, COLORREF color) {
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);
    MoveToEx(hdc, xmin, ymin, NULL);
    LineTo(hdc, xmax, ymin);
    LineTo(hdc, xmax, ymax);
    LineTo(hdc, xmin, ymax);
    LineTo(hdc, xmin, ymin);
    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

// Cohen–Sutherland outcodes
enum {
    INSIDE_CS = 0, // 0000
    LEFT_CS   = 1, // 0001
    RIGHT_CS  = 2, // 0010
    BOTTOM_CS = 4, // 0100
    TOP_CS    = 8  // 1000
};

static int computeOutCode(int x, int y, int xmin, int ymin, int xmax, int ymax) {
    int code = INSIDE_CS;
    if (x < xmin)       code |= LEFT_CS;
    else if (x > xmax)  code |= RIGHT_CS;
    if (y < ymin)       code |= BOTTOM_CS;
    else if (y > ymax)  code |= TOP_CS;
    return code;
}

static bool cohenSutherlandClip(int x1, int y1, int x2, int y2,
                                int xmin, int ymin, int xmax, int ymax,
                                int &cx1, int &cy1, int &cx2, int &cy2)
{
    int out1 = computeOutCode(x1, y1, xmin, ymin, xmax, ymax);
    int out2 = computeOutCode(x2, y2, xmin, ymin, xmax, ymax);
    bool accept = false;
    double dx = x2 - x1, dy = y2 - y1;

    while (true) {
        if ((out1 | out2) == 0) {
            // both inside
            cx1 = x1; cy1 = y1;
            cx2 = x2; cy2 = y2;
            accept = true;
            break;
        }
        if (out1 & out2) {
            // both outside same region
            break;
        }
        int outcodeOut = out1 ? out1 : out2;
        double x, y;

        if (outcodeOut & TOP_CS) {
            x = x1 + dx * (ymax - y1) / dy;  y = ymax;
        }
        else if (outcodeOut & BOTTOM_CS) {
            x = x1 + dx * (ymin - y1) / dy;  y = ymin;
        }
        else if (outcodeOut & RIGHT_CS) {
            y = y1 + dy * (xmax - x1) / dx;  x = xmax;
        }
        else { // LEFT_CS
            y = y1 + dy * (xmin - x1) / dx;  x = xmin;
        }

        if (outcodeOut == out1) {
            x1 = (int)round(x); y1 = (int)round(y);
            out1 = computeOutCode(x1, y1, xmin, ymin, xmax, ymax);
        } else {
            x2 = (int)round(x); y2 = (int)round(y);
            out2 = computeOutCode(x2, y2, xmin, ymin, xmax, ymax);
        }
    }
    return accept;
}

// Sutherland–Hodgman polygon clipping
static POINT intersectPolyEdge(const POINT &P1, const POINT &P2, int edge,
                               int xmin, int ymin, int xmax, int ymax)
{
    double x1 = P1.x, y1 = P1.y;
    double x2 = P2.x, y2 = P2.y;
    double dx = x2 - x1, dy = y2 - y1;
    double x, y;

    switch (edge) {
        case 0: // LEFT
            x = xmin;
            y = y1 + dy * (xmin - x1) / dx;
            break;
        case 1: // RIGHT
            x = xmax;
            y = y1 + dy * (xmax - x1) / dx;
            break;
        case 2: // BOTTOM
            y = ymin;
            x = x1 + dx * (ymin - y1) / dy;
            break;
        case 3: // TOP
            y = ymax;
            x = x1 + dx * (ymax - y1) / dy;
            break;
        default:
            x = x1; y = y1;
    }
    return { (int)round(x), (int)round(y) };
}

static bool insidePolyEdge(const POINT &P, int edge,
                           int xmin, int ymin, int xmax, int ymax)
{
    switch (edge) {
        case 0: return P.x >= xmin;  // LEFT
        case 1: return P.x <= xmax;  // RIGHT
        case 2: return P.y >= ymin;  // BOTTOM
        case 3: return P.y <= ymax;  // TOP
    }
    return false;
}

static void sutherlandHodgmanPolygonClip(const vector<POINT> &inPoly,
                                         vector<POINT> &outPoly,
                                         int xmin, int ymin,
                                         int xmax, int ymax)
{
    vector<POINT> tempPoly = inPoly;

    for (int edge = 0; edge < 4; edge++) {
        outPoly.clear();
        if (tempPoly.empty()) break;
        POINT S = tempPoly.back();

        for (const POINT &E : tempPoly) {
            bool insideE = insidePolyEdge(E, edge, xmin, ymin, xmax, ymax);
            bool insideS = insidePolyEdge(S, edge, xmin, ymin, xmax, ymax);

            if (insideE) {
                if (!insideS) {
                    POINT I = intersectPolyEdge(S, E, edge, xmin, ymin, xmax, ymax);
                    outPoly.push_back(I);
                }
                outPoly.push_back(E);
            }
            else if (insideS) {
                POINT I = intersectPolyEdge(S, E, edge, xmin, ymin, xmax, ymax);
                outPoly.push_back(I);
            }
            S = E;
        }
        tempPoly = outPoly;
    }
}

//-----------------------------------------------------------------------------
// Requirement 20: Square Clipping [Point, Line]
//-----------------------------------------------------------------------------

// Reuse computeOutCode and cohenSutherlandClip by treating square as rectangle

//-----------------------------------------------------------------------------
// Bonus: Circular Clipping [Point, Line]
//-----------------------------------------------------------------------------

static void pointClipCircle(HDC hdc, int x, int y,
                            int radius, int centerX, int centerY, COLORREF color)
{
    double dx2 = pow(centerX - x, 2);
    double dy2 = pow(centerY - y, 2);
    double dist = sqrt(dx2 + dy2);
    if (dist <= radius) {
        SetPixel(hdc, x, y, RGB(0, 0, 255)); // inside = blue
    } else {
        SetPixel(hdc, x, y, RGB(255, 0, 0)); // outside = red
    }
}

static void circleClipLine(HDC hdc, int x1, int y1, int x2, int y2,
                           int cx, int cy, int r, COLORREF color)
{
    // Parametrize line: P(t) = (x1, y1) + t * (dx, dy), t in [0,1]
    double dx = x2 - x1;
    double dy = y2 - y1;
    double fx = x1 - cx;
    double fy = y1 - cy;

    double A = dx * dx + dy * dy;
    double B = 2.0 * (dx * fx + dy * fy);
    double C = fx * fx + fy * fy - (double)r * r;
    double disc = B * B - 4.0 * A * C;

    bool inside1 = (fx * fx + fy * fy <= (double)r * r);
    double t0 = 0.0, t1 = 1.0;

    if (disc < 0.0) {
        // no intersections
        if (inside1) {
            // entire segment inside
            lineBresenham(hdc, x1, y1, x2, y2, color);
            shapes.push_back({ "ClipCircleLine", x1, y1, x2, y2, cx, cy, r, color });
        }
        // else: entirely outside, draw nothing
        return;
    }

    double sqrt_disc = sqrt(disc);
    double tA = (-B + sqrt_disc) / (2.0 * A);
    double tB = (-B - sqrt_disc) / (2.0 * A);
    t0 = min(tA, tB);
    t1 = max(tA, tB);

    vector<pair<double, POINT>> pts;
    if (inside1) {
        pts.push_back({ 0.0,      { x1, y1 } });
    }
    if (t0 > 0.0 && t0 < 1.0) {
        int xi = (int)round(x1 + t0 * dx);
        int yi = (int)round(y1 + t0 * dy);
        pts.push_back({ t0, { xi, yi } });
    }
    if (t1 > 0.0 && t1 < 1.0) {
        int xi = (int)round(x1 + t1 * dx);
        int yi = (int)round(y1 + t1 * dy);
        pts.push_back({ t1, { xi, yi } });
    }
    double fx2 = x2 - cx;
    double fy2 = y2 - cy;
    bool inside2 = (fx2 * fx2 + fy2 * fy2 <= (double)r * r);
    if (!inside1 && inside2) {
        pts.push_back({ 1.0, { x2, y2 } });
    }

    if (pts.size() < 2) {
        // no valid inside segment
        return;
    }

    sort(pts.begin(), pts.end(),
         [](auto &a, auto &b) { return a.first < b.first; });

    // if three points (one endpoint inside), take middle two
    int startIdx = (pts.size() == 3 ? 1 : 0);
    POINT P0 = pts[startIdx].second;
    POINT P1 = pts[startIdx + 1].second;
    lineBresenham(hdc, P0.x, P0.y, P1.x, P1.y, color);
    shapes.push_back({ "ClipCircleLine", P0.x, P0.y, P1.x, P1.y, cx, cy, r, color });
}

//-----------------------------------------------------------------------------
// File I/O: Save & Load
//-----------------------------------------------------------------------------

static void PerformSave(const char* fname) {
    ofstream out(fname);
    if (!out.is_open()) return;
    for (auto &s : shapes) {
        out << s.type << " "
            << s.x1 << " " << s.y1 << " "
            << s.x2 << " " << s.y2 << " "
            << s.a  << " " << s.b  << " "
            << s.c  << " "
            << s.color;
        if (!s.verts.empty()) {
            out << " " << s.verts.size();
            for (auto &p : s.verts) {
                out << " " << p.x << " " << p.y;
            }
        }
        out << "\n";
    }
    out.close();
}

static void PerformLoad(const char* fname) {
    ifstream in(fname);
    if (!in.is_open()) return;
    shapes.clear();
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        Shape s;
        iss >> s.type
            >> s.x1 >> s.y1
            >> s.x2 >> s.y2
            >> s.a  >> s.b
            >> s.c
            >> s.color;
        if (s.type == "ConvexPoly" ||
            s.type == "NonConvexPoly" ||
            s.type == "CardinalSpline" ||
            s.type == "ClipRectPoly" ||
            s.type == "ClipCirclePoly")
        {
            size_t n;
            iss >> n;
            s.verts.resize(n);
            for (size_t i = 0; i < n; i++) {
                iss >> s.verts[i].x >> s.verts[i].y;
            }
        }
        shapes.push_back(s);
    }
    in.close();
    InvalidateRect(GetActiveWindow(), NULL, TRUE);
}

//-----------------------------------------------------------------------------
// Menu Creation (grouped by requirement categories)
//-----------------------------------------------------------------------------

static void CreateMenus(HWND hwnd) {
    hMenu = CreateMenu();

    // ------------- File & Clear (Req 6, 7, 8) -------------
    HMENU hFile = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFile, "File");
    AppendMenu(hFile, MF_STRING, ID_SAVE,  "Save Shapes");   // Req 7
    AppendMenu(hFile, MF_STRING, ID_LOAD,  "Load Shapes");   // Req 8
    AppendMenu(hFile, MF_SEPARATOR, 0, NULL);
    AppendMenu(hFile, MF_STRING, ID_CLEAR, "Clear Screen");  // Req 6

    // ------------- Color Selection (Req 5) -------------
    HMENU hColor = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hColor, "Color");
    AppendMenu(hColor, MF_STRING, ID_COLOR_RED,   "Red");
    AppendMenu(hColor, MF_STRING, ID_COLOR_BLACK, "Black");
    AppendMenu(hColor, MF_STRING, ID_COLOR_BLUE,  "Blue");
    AppendMenu(hColor, MF_STRING, ID_COLOR_WHITE, "White");
    AppendMenu(hColor, MF_STRING, ID_COLOR_PINK,  "Pink");

    // ------------- Drawing Algorithms (Req 9, 10, 18) -------------
    // Line algorithms (Req 9)
    HMENU hLine = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hLine, "Line Algorithms");
    AppendMenu(hLine, MF_STRING, ID_LINE_DDA,   "DDA");       // Req 9
    AppendMenu(hLine, MF_STRING, ID_LINE_MID,   "Midpoint");  // Req 9
    AppendMenu(hLine, MF_STRING, ID_LINE_PARAM, "Parametric"); // Req 9

    // Circle algorithms (Req 10)
    HMENU hCircle = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCircle, "Circle Algorithms");
    AppendMenu(hCircle, MF_STRING, ID_CIRCLE_DIRECT,    "Direct");           // Req 10
    AppendMenu(hCircle, MF_STRING, ID_CIRCLE_POLAR,     "Polar");            // Req 10
    AppendMenu(hCircle, MF_STRING, ID_CIRCLE_ITERATIVE, "Iterative Polar");  // Req 10
    AppendMenu(hCircle, MF_STRING, ID_CIRCLE_MIDPOINT,  "Midpoint");         // Req 10
    AppendMenu(hCircle, MF_STRING, ID_CIRCLE_FAST,      "Fast Bresenham");   // Req 10

    // Ellipse algorithms (Req 18)
    HMENU hEllipse = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hEllipse, "Ellipse Algorithms");
    AppendMenu(hEllipse, MF_STRING, ID_ELLIPSE_DIRECT,    "Direct");   // Req 18
    AppendMenu(hEllipse, MF_STRING, ID_ELLIPSE_POLAR,     "Polar");    // Req 18
    AppendMenu(hEllipse, MF_STRING, ID_ELLIPSE_MIDPOINT,  "Midpoint"); // Req 18

    // ------------- Filling Algorithms (Req 11, 12, 13, 14, 15, 16, 17) -------------
    // Circle Filling (lines) (Req 11)
    HMENU hFillCircleLines = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFillCircleLines, "Fill Circle with Lines");
    AppendMenu(hFillCircleLines, MF_STRING, ID_FILL_CIRCLE_Q1, "Quadrant 1");
    AppendMenu(hFillCircleLines, MF_STRING, ID_FILL_CIRCLE_Q2, "Quadrant 2");
    AppendMenu(hFillCircleLines, MF_STRING, ID_FILL_CIRCLE_Q3, "Quadrant 3");
    AppendMenu(hFillCircleLines, MF_STRING, ID_FILL_CIRCLE_Q4, "Quadrant 4");

    // Circle Filling (circles) (Req 12)
    HMENU hFillCircleCircles = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFillCircleCircles, "Fill Circle with Circles");
    AppendMenu(hFillCircleCircles, MF_STRING, ID_FILL_LINE_Q1, "Quadrant 1");
    AppendMenu(hFillCircleCircles, MF_STRING, ID_FILL_LINE_Q2, "Quadrant 2");
    AppendMenu(hFillCircleCircles, MF_STRING, ID_FILL_LINE_Q3, "Quadrant 3");
    AppendMenu(hFillCircleCircles, MF_STRING, ID_FILL_LINE_Q4, "Quadrant 4");

    // Square Filling with Hermite (Req 13)
    // Changed to a direct menu item instead of a submenu
    AppendMenu(hMenu, MF_STRING, ID_HERMITE_SQUARE, "Fill Square (Hermite)");

    // Rectangle Filling with Bezier (Req 14)
    HMENU hBezier = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hBezier, "Fill Rectangle (Bezier)");
    AppendMenu(hBezier, MF_STRING, ID_FILL_BEZIERRECT, "Horizontal"); // Req 14

    // Polygon Filling (Req 15)
    HMENU hPolyFill = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hPolyFill, "Fill Polygon");
    AppendMenu(hPolyFill, MF_STRING, ID_FILL_CONVEX_POLY,    "Convex");
    AppendMenu(hPolyFill, MF_STRING, ID_FILL_NONCONVEX_POLY, "Non-Convex");

    // Flood Fill (Req 16)
    HMENU hFlood = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFlood, "Flood Fill");
    AppendMenu(hFlood, MF_STRING, ID_FLOOD_REC,   "Recursive");
    AppendMenu(hFlood, MF_STRING, ID_FLOOD_STACK, "Non-Recursive");

    // Cardinal Spline (Req 17)
    HMENU hSpline = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hSpline, "Curves");
    AppendMenu(hSpline, MF_STRING, ID_CARDINAL_SPLINE, "Cardinal Spline");

    // ------------- Clipping Algorithms (Req 19, 20, Bonus) -------------
    // Rectangular clipping (Req 19)
    HMENU hClipRect = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClipRect, "Rectangular Clipping");
    AppendMenu(hClipRect, MF_STRING, ID_CLIP_RECT_WINDOW, "Define Window");
    AppendMenu(hClipRect, MF_STRING, ID_CLIP_RECT_POINT,  "Clip Point");
    AppendMenu(hClipRect, MF_STRING, ID_CLIP_RECT_LINE,   "Clip Line");
    AppendMenu(hClipRect, MF_STRING, ID_CLIP_RECT_POLY,   "Clip Polygon");

    // Square clipping (Req 20)
    HMENU hClipSquare = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClipSquare, "Square Clipping");
    AppendMenu(hClipSquare, MF_STRING, ID_CLIP_SQUARE_WINDOW, "Define Window");
    AppendMenu(hClipSquare, MF_STRING, ID_CLIP_SQUARE_POINT,  "Clip Point");
    AppendMenu(hClipSquare, MF_STRING, ID_CLIP_SQUARE_LINE,   "Clip Line");

    // Circular clipping (Bonus)
    HMENU hClipCircle = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClipCircle, "Circular Clipping");
    AppendMenu(hClipCircle, MF_STRING, ID_CLIP_CIRCLE_WINDOW, "Define Window");
    AppendMenu(hClipCircle, MF_STRING, ID_CLIP_CIRCLE_POINT,  "Clip Point");
    AppendMenu(hClipCircle, MF_STRING, ID_CLIP_CIRCLE_LINE,   "Clip Line");

    SetMenu(hwnd, hMenu);
}

//-----------------------------------------------------------------------------
// Window Procedure & Message Handling
//-----------------------------------------------------------------------------

LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
    static int currentMode = 0;
    static COLORREF currentColor = RGB(0, 0, 0);

    HDC hdc = GetDC(hwnd);

    switch (message) {
    case WM_COMMAND:
        switch (wParam) {
            // File operations
            case ID_SAVE:
                PerformSave("shapes_data.txt");
                cout << "Saved shapes to file\n";
                break;
            case ID_LOAD:
                PerformLoad("shapes_data.txt");
                cout << "Loaded shapes from file\n";
                break;
            case ID_CLEAR:
                shapes.clear();
                currentVerts.clear();
                InvalidateRect(hwnd, NULL, TRUE);
                cout << "Cleared canvas\n";
                break;

            // Color selection
            case ID_COLOR_RED:    currentColor = RGB(255, 0, 0);    cout << "Color: Red\n";    break;
            case ID_COLOR_BLACK:  currentColor = RGB(0, 0, 0);      cout << "Color: Black\n";  break;
            case ID_COLOR_BLUE:   currentColor = RGB(0, 0, 255);    cout << "Color: Blue\n";   break;
            case ID_COLOR_WHITE:  currentColor = RGB(255, 255, 255);cout << "Color: White\n";  break;
            case ID_COLOR_PINK:   currentColor = RGB(255, 63, 127); cout << "Color: Pink\n";   break;

            // Line algorithms
            case ID_LINE_DDA:
                currentMode = ID_LINE_DDA;
                cout << "Mode: DDA Line\n";
                break;
            case ID_LINE_MID:
                currentMode = ID_LINE_MID;
                cout << "Mode: Midpoint Line\n";
                break;
            case ID_LINE_PARAM:
                currentMode = ID_LINE_PARAM;
                cout << "Mode: Parametric Line\n";
                break;

            // Circle algorithms
            case ID_CIRCLE_DIRECT:
                currentMode = ID_CIRCLE_DIRECT;
                cout << "Mode: Circle Direct\n";
                break;
            case ID_CIRCLE_POLAR:
                currentMode = ID_CIRCLE_POLAR;
                cout << "Mode: Circle Polar\n";
                break;
            case ID_CIRCLE_ITERATIVE:
                currentMode = ID_CIRCLE_ITERATIVE;
                cout << "Mode: Circle Iterative Polar\n";
                break;
            case ID_CIRCLE_MIDPOINT:
                currentMode = ID_CIRCLE_MIDPOINT;
                cout << "Mode: Circle Midpoint\n";
                break;
            case ID_CIRCLE_FAST:
                currentMode = ID_CIRCLE_FAST;
                cout << "Mode: Circle Fast Bresenham\n";
                break;

            // Filling circle with lines
            case ID_FILL_CIRCLE_Q1:
            case ID_FILL_CIRCLE_Q2:
            case ID_FILL_CIRCLE_Q3:
            case ID_FILL_CIRCLE_Q4: {
                int quad = (int)wParam - ID_FILL_CIRCLE_Q1 + 1; // 1..4
                bool found = false;
                int cx = 0, cy = 0, radius = 0;
                for (int i = (int)shapes.size() - 1; i >= 0; i--) {
                    auto &s = shapes[i];
                    if (s.type == "CircleDirect" ||
                        s.type == "CirclePolar"  ||
                        s.type == "CircleIterative" ||
                        s.type == "CircleMidpoint" ||
                        s.type == "CircleFast")
                    {
                        found = true;
                        cx = s.x1;
                        cy = s.y1;
                        radius = s.a;
                        break;
                    }
                }
                if (found) {
                    fillWithLines(hdc, cx, cy, 0, 0, radius, currentColor, quad);
                    shapes.push_back({ "FillCircleLines", cx, cy, 0, 0, radius, quad, 0, currentColor });
                    cout << "Filled last circle with lines in quadrant " << quad << "\n";
                } else {
                    cout << "No circle found to fill!\n";
                }
                InvalidateRect(hwnd, NULL, TRUE);
                break;
            }

            // Filling circle with circles
            case ID_FILL_LINE_Q1:
            case ID_FILL_LINE_Q2:
            case ID_FILL_LINE_Q3:
            case ID_FILL_LINE_Q4: {
                int quad = (int)wParam - ID_FILL_LINE_Q1 + 1; // 1..4
                bool found = false;
                int cx = 0, cy = 0, radius = 0;
                for (int i = (int)shapes.size() - 1; i >= 0; i--) {
                    auto &s = shapes[i];
                    if (s.type == "CircleDirect" ||
                        s.type == "CirclePolar"  ||
                        s.type == "CircleIterative" ||
                        s.type == "CircleMidpoint" ||
                        s.type == "CircleFast")
                    {
                        found = true;
                        cx = s.x1;
                        cy = s.y1;
                        radius = s.a;
                        break;
                    }
                }
                if (found) {
                    fillWithCircles(hdc, cx, cy, 0, 0, radius, currentColor, quad);
                    shapes.push_back({ "FillCircleCircles", cx, cy, 0, 0, radius, quad, 0, currentColor });
                    cout << "Filled last circle with circles in quadrant " << quad << "\n";
                } else {
                    cout << "No circle found to fill!\n";
                }
                InvalidateRect(hwnd, NULL, TRUE);
                break;
            }

            // Fill Square (Hermite)
            case ID_HERMITE_SQUARE:
                currentMode = ID_HERMITE_SQUARE;
                cout << "Mode: Fill Square with Hermite\n";
                break;

            // Bezier rectangle fill
            case ID_FILL_BEZIERRECT:
                currentMode = ID_FILL_BEZIERRECT;
                cout << "Mode: Fill Rectangle with Bezier\n";
                break;

            // Polygon fills
            case ID_FILL_CONVEX_POLY:
                currentMode = ID_FILL_CONVEX_POLY;
                currentVerts.clear();
                cout << "Mode: Fill Convex Polygon\n";
                break;
            case ID_FILL_NONCONVEX_POLY:
                currentMode = ID_FILL_NONCONVEX_POLY;
                currentVerts.clear();
                cout << "Mode: Fill Non-Convex Polygon\n";
                break;

            // Cardinal spline
            case ID_CARDINAL_SPLINE:
                currentMode = ID_CARDINAL_SPLINE;
                currentVerts.clear();
                cout << "Mode: Cardinal Spline\n";
                break;

            // Flood fill
            case ID_FLOOD_REC:
                currentMode = ID_FLOOD_REC;
                cout << "Mode: Flood Fill (Recursive)\n";
                break;
            case ID_FLOOD_STACK:
                currentMode = ID_FLOOD_STACK;
                cout << "Mode: Flood Fill (Non-Recursive)\n";
                break;

            // Ellipse algorithms
            case ID_ELLIPSE_DIRECT:
                currentMode = ID_ELLIPSE_DIRECT;
                cout << "Mode: Ellipse Direct\n";
                break;
            case ID_ELLIPSE_POLAR:
                currentMode = ID_ELLIPSE_POLAR;
                cout << "Mode: Ellipse Polar\n";
                break;
            case ID_ELLIPSE_MIDPOINT:
                currentMode = ID_ELLIPSE_MIDPOINT;
                cout << "Mode: Ellipse Midpoint\n";
                break;

            // Rectangular clipping
            case ID_CLIP_RECT_WINDOW:
                currentMode = ID_CLIP_RECT_WINDOW;
                cout << "Mode: Define Rectangular Clip Window\n";
                break;
            case ID_CLIP_RECT_POINT:
                currentMode = ID_CLIP_RECT_POINT;
                cout << "Mode: Clip Point in Rectangular Window\n";
                break;
            case ID_CLIP_RECT_LINE:
                currentMode = ID_CLIP_RECT_LINE;
                cout << "Mode: Clip Line in Rectangular Window\n";
                break;
            case ID_CLIP_RECT_POLY:
                currentMode = ID_CLIP_RECT_POLY;
                currentVerts.clear();
                cout << "Mode: Clip Polygon in Rectangular Window\n";
                break;

            // Square clipping
            case ID_CLIP_SQUARE_WINDOW:
                currentMode = ID_CLIP_SQUARE_WINDOW;
                cout << "Mode: Define Square Clip Window\n";
                break;
            case ID_CLIP_SQUARE_POINT:
                currentMode = ID_CLIP_SQUARE_POINT;
                cout << "Mode: Clip Point in Square Window\n";
                break;
            case ID_CLIP_SQUARE_LINE:
                currentMode = ID_CLIP_SQUARE_LINE;
                cout << "Mode: Clip Line in Square Window\n";
                break;

            // Circular clipping (bonus)
            case ID_CLIP_CIRCLE_WINDOW:
                currentMode = ID_CLIP_CIRCLE_WINDOW;
                cout << "Mode: Define Circular Clip Window\n";
                break;
            case ID_CLIP_CIRCLE_LINE:
                currentMode = ID_CLIP_CIRCLE_LINE;
                cout << "Mode: Clip Line in Circular Window\n";
                break;
            case ID_CLIP_CIRCLE_POINT:
                currentMode = ID_CLIP_CIRCLE_POINT;
                cout << "Mode: Clip Point in Circular Window\n";
                break;

            default:
                break;
        }
        break;

    case WM_LBUTTONDOWN:
        xc = LOWORD(lParam);
        yc = HIWORD(lParam);

        if (currentMode == ID_CLIP_RECT_WINDOW) {
            rectXmin = xc;
            rectYmin = yc;
            cout << "Rectangular clip window first corner: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_CLIP_SQUARE_WINDOW) {
            sqXmin = xc;
            sqYmin = yc;
            cout << "Square clip window first corner: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_CLIP_CIRCLE_WINDOW) {
            windowX = xc;
            windowY = yc;
            cout << "Circular clip center set to: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_FILL_CONVEX_POLY) {
            currentVerts.push_back({ xc, yc });
            SetPixel(hdc, xc, yc, currentColor);
            cout << "Added convex-polygon vertex: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_FILL_NONCONVEX_POLY) {
            currentVerts.push_back({ xc, yc });
            SetPixel(hdc, xc, yc, currentColor);
            cout << "Added non-convex-polygon vertex: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_CARDINAL_SPLINE) {
            currentVerts.push_back({ xc, yc });
            SetPixel(hdc, xc, yc, currentColor);
            cout << "Added spline control point: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_HERMITE_SQUARE) {
            cout << "Set Hermite square top-left: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_FLOOD_REC || currentMode == ID_FLOOD_STACK) {
            cout << "Flood fill seed: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_FILL_BEZIERRECT) {
            cout << "Set Bezier rectangle top-left: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_CLIP_RECT_POINT) {
            cout << "Selected rectangular clip point: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_CLIP_SQUARE_POINT) {
            cout << "Selected square clip point: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_CLIP_CIRCLE_POINT) {
            cout << "Selected circular clip point: (" << xc << "," << yc << ")\n";
        }
        else if (currentMode == ID_CLIP_RECT_POLY) {
            currentVerts.push_back({ xc, yc });
            SetPixel(hdc, xc, yc, currentColor);
            cout << "Added clip-polygon vertex: (" << xc << "," << yc << ")\n";
        }
        break;

    case WM_RBUTTONDOWN:
        xe = LOWORD(lParam);
        ye = HIWORD(lParam);

        switch (currentMode) {
            // Draw lines
            case ID_LINE_DDA: {
                DDALine(hdc, xc, yc, xe, ye, currentColor);
                shapes.push_back({ "DDA", xc, yc, xe, ye, 0, 0, 0, currentColor });
                cout << "Drew DDA line: (" << xc << "," << yc << ") -> (" << xe << "," << ye << ")\n";
                break;
            }
            case ID_LINE_MID: {
                lineBresenham(hdc, xc, yc, xe, ye, currentColor);
                shapes.push_back({ "MidpointLine", xc, yc, xe, ye, 0, 0, 0, currentColor });
                cout << "Drew Midpoint line: (" << xc << "," << yc << ") -> (" << xe << "," << ye << ")\n";
                break;
            }
            case ID_LINE_PARAM: {
                parametricLine(hdc, xc, yc, xe, ye, 0, 0, 0, currentColor);
                shapes.push_back({ "ParametricLine", xc, yc, xe, ye, 0, 0, 0, currentColor });
                cout << "Drew Parametric line: (" << xc << "," << yc << ") -> (" << xe << "," << ye << ")\n";
                break;
            }

            // Draw circles
            case ID_CIRCLE_DIRECT: {
                R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
                circleDirectMethod(hdc, xc, yc, R, currentColor);
                shapes.push_back({ "CircleDirect", xc, yc, 0, 0, R, 0, 0, currentColor });
                cout << "Drew Direct circle: center(" << xc << "," << yc << "), r=" << R << "\n";
                break;
            }
            case ID_CIRCLE_POLAR: {
                R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
                circlePolar(hdc, xc, yc, R, currentColor);
                shapes.push_back({ "CirclePolar", xc, yc, 0, 0, R, 0, 0, currentColor });
                cout << "Drew Polar circle: center(" << xc << "," << yc << "), r=" << R << "\n";
                break;
            }
            case ID_CIRCLE_ITERATIVE: {
                R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
                circleIterative(hdc, xc, yc, R, currentColor);
                shapes.push_back({ "CircleIterative", xc, yc, 0, 0, R, 0, 0, currentColor });
                cout << "Drew Iterative Polar circle: center(" << xc << "," << yc << "), r=" << R << "\n";
                break;
            }
            case ID_CIRCLE_MIDPOINT: {
                R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
                midpointCircle(hdc, xc, yc, R, currentColor);
                shapes.push_back({ "CircleMidpoint", xc, yc, 0, 0, R, 0, 0, currentColor });
                cout << "Drew Midpoint circle: center(" << xc << "," << yc << "), r=" << R << "\n";
                break;
            }
            case ID_CIRCLE_FAST: {
                R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
                circleFastBresenham(hdc, xc, yc, R, currentColor);
                shapes.push_back({ "CircleFast", xc, yc, 0, 0, R, 0, 0, currentColor });
                cout << "Drew Fast Bresenham circle: center(" << xc << "," << yc << "), r=" << R << "\n";
                break;
            }

            // Fill Square (Hermite)
            case ID_HERMITE_SQUARE: {
                POINT tl = { xc, yc };
                int dx = abs(xe - xc), dy = abs(ye - yc);
                int size = min(dx, dy);
                FillHermiteSquare(hdc, tl, size, currentColor);
                shapes.push_back({ "HermiteSquare", xc, yc, 0, 0, size, 0, 0, currentColor });
                cout << "Filled square with Hermite at (" << xc << "," << yc << "), size=" << size << "\n";
                break;
            }

            // Bezier rectangle fill
            case ID_FILL_BEZIERRECT: {
                int left   = min(xc, xe);
                int right  = max(xc, xe);
                int top    = min(yc, ye);
                int bottom = max(yc, ye);
                POINT tl = { left, top };
                POINT br = { right, bottom };
                FillBezierRectangle(hdc, tl, br, currentColor);
                shapes.push_back({ "FillBezierRect", left, top, right, bottom, 0, 0, 0, currentColor });
                cout << "Filled rectangle with Bezier: (" << left << "," << top << ") -> (" << right << "," << bottom << ")\n";
                break;
            }

            // Polygon fills
            case ID_FILL_CONVEX_POLY: {
                if (currentVerts.size() >= 3) {
                    fillConvexPolygon(hdc, currentVerts, currentColor);
                    Shape poly = { "ConvexPoly", 0, 0, 0, 0, 0, 0, 0, currentColor, currentVerts };
                    shapes.push_back(poly);
                    cout << "Filled convex polygon with " << currentVerts.size() << " vertices\n";
                } else {
                    cout << "Not enough vertices for convex polygon\n";
                }
                currentVerts.clear();
                break;
            }
            case ID_FILL_NONCONVEX_POLY: {
                if (currentVerts.size() >= 3) {
                    fillNonConvexPolygon(hdc, currentVerts, currentColor);
                    Shape poly = { "NonConvexPoly", 0, 0, 0, 0, 0, 0, 0, currentColor, currentVerts };
                    shapes.push_back(poly);
                    cout << "Filled non-convex polygon with " << currentVerts.size() << " vertices\n";
                } else {
                    cout << "Not enough vertices for non-convex polygon\n";
                }
                currentVerts.clear();
                break;
            }

            // Cardinal spline
            case ID_CARDINAL_SPLINE: {
                if (currentVerts.size() >= 3) {
                    drawCardinalSpline(hdc, currentVerts, currentColor);
                    Shape spline = { "CardinalSpline", 0, 0, 0, 0, 0, 0, 0, currentColor, currentVerts };
                    shapes.push_back(spline);
                    cout << "Drew Cardinal spline with " << currentVerts.size() << " control points\n";
                } else {
                    cout << "Not enough control points for spline\n";
                }
                currentVerts.clear();
                break;
            }

            // Flood fill
            case ID_FLOOD_REC: {
                COLORREF tar = GetPixel(hdc, xc, yc);
                if (tar != currentColor) {
                    floodFillRec(hdc, xc, yc, tar, currentColor);
                    shapes.push_back({ "FloodRec", xc, yc, 0, 0, 0, 0, 0, currentColor });
                    cout << "Performed recursive flood-fill at (" << xc << "," << yc << ")\n";
                } else {
                    cout << "Seed color equals fill color\n";
                }
                break;
            }
            case ID_FLOOD_STACK: {
                COLORREF tar = GetPixel(hdc, xc, yc);
                if (tar != currentColor) {
                    floodFillStack(hdc, xc, yc, tar, currentColor);
                    shapes.push_back({ "FloodStack", xc, yc, 0, 0, 0, 0, 0, currentColor });
                    cout << "Performed non-recursive flood-fill at (" << xc << "," << yc << ")\n";
                } else {
                    cout << "Seed color equals fill color\n";
                }
                break;
            }

            // Ellipse algorithms
            case ID_ELLIPSE_DIRECT: {
                int a = abs(xe - xc), b = abs(ye - yc);
                ellipseDirect(hdc, xc, yc, a, b, currentColor);
                shapes.push_back({ "EllipseDirect", xc, yc, 0, 0, a, b, 0, currentColor });
                cout << "Drew Direct ellipse: center(" << xc << "," << yc << "), a=" << a << ", b=" << b << "\n";
                break;
            }
            case ID_ELLIPSE_POLAR: {
                int a = abs(xe - xc), b = abs(ye - yc);
                ellipsePolar(hdc, xc, yc, a, b, currentColor);
                shapes.push_back({ "EllipsePolar", xc, yc, 0, 0, a, b, 0, currentColor });
                cout << "Drew Polar ellipse: center(" << xc << "," << yc << "), a=" << a << ", b=" << b << "\n";
                break;
            }
            case ID_ELLIPSE_MIDPOINT: {
                int a = abs(xe - xc), b = abs(ye - yc);
                ellipseMidpoint(hdc, xc, yc, a, b, currentColor);
                shapes.push_back({ "EllipseMidpoint", xc, yc, 0, 0, a, b, 0, currentColor });
                cout << "Drew Midpoint ellipse: center(" << xc << "," << yc << "), a=" << a << ", b=" << b << "\n";
                break;
            }

            // Rectangular clipping
            case ID_CLIP_RECT_WINDOW: {
                rectXmin = xc;
                rectYmin = yc;
                rectXmax = xe;
                rectYmax = ye;
                if (rectXmax < rectXmin) swap(rectXmin, rectXmax);
                if (rectYmax < rectYmin) swap(rectYmin, rectYmax);
                drawRectangleWindow(hdc, rectXmin, rectYmin, rectXmax, rectYmax, currentColor);
                shapes.push_back({ "ClipRectWindow", rectXmin, rectYmin, rectXmax, rectYmax, 0, 0, 0, currentColor });
                cout << "Rectangular clipping window: (" << rectXmin << "," << rectYmin << ") -> (" << rectXmax << "," << rectYmax << ")\n";
                break;
            }
            case ID_CLIP_RECT_POINT: {
                bool inside = (xc >= rectXmin && xc <= rectXmax && yc >= rectYmin && yc <= rectYmax);
                SetPixel(hdc, xc, yc, inside ? RGB(0,0,255) : RGB(255,0,0));
                shapes.push_back({ "ClipRectPoint", xc, yc, 0, 0, 0, 0, 0, inside ? RGB(0,0,255) : RGB(255,0,0) });
                cout << "Rectangular clip point (" << xc << "," << yc << ") is " << (inside ? "inside\n" : "outside\n");
                break;
            }
            case ID_CLIP_RECT_LINE: {
                int cx1, cy1, cx2, cy2;
                if (cohenSutherlandClip(xc, yc, xe, ye, rectXmin, rectYmin, rectXmax, rectYmax, cx1, cy1, cx2, cy2)) {
                    lineBresenham(hdc, cx1, cy1, cx2, cy2, currentColor);
                    shapes.push_back({ "ClipRectLine", cx1, cy1, cx2, cy2, rectXmin, rectYmin, 0, currentColor });
                    cout << "Clipped line in rectangular window: (" << cx1 << "," << cy1 << ") -> (" << cx2 << "," << cy2 << ")\n";
                } else {
                    cout << "Line entirely outside rectangular window\n";
                }
                break;
            }
            case ID_CLIP_RECT_POLY: {
                if (currentVerts.size() >= 3) {
                    vector<POINT> clipped;
                    sutherlandHodgmanPolygonClip(currentVerts, clipped, rectXmin, rectYmin, rectXmax, rectYmax);
                    if (!clipped.empty()) {
                        HPEN hPenC = CreatePen(PS_SOLID, 1, currentColor);
                        HGDIOBJ oldPenC = SelectObject(hdc, hPenC);
                        MoveToEx(hdc, clipped[0].x, clipped[0].y, NULL);
                        for (size_t i = 1; i < clipped.size(); i++) {
                            LineTo(hdc, clipped[i].x, clipped[i].y);
                        }
                        LineTo(hdc, clipped[0].x, clipped[0].y);
                        SelectObject(hdc, oldPenC);
                        DeleteObject(hPenC);

                        Shape poly = { "ClipRectPoly", 0, 0, 0, 0, 0, 0, 0, currentColor, clipped };
                        shapes.push_back(poly);
                        cout << "Clipped polygon in rectangular window\n";
                    } else {
                        cout << "Clipped polygon empty (fully outside)\n";
                    }
                } else {
                    cout << "Not enough vertices for polygon clipping\n";
                }
                currentVerts.clear();
                break;
            }

            // Square clipping
            case ID_CLIP_SQUARE_WINDOW: {
                sqXmin = xc;
                sqYmin = yc;
                int dx = xe - xc;
                int dy = ye - yc;
                sqSide = max(abs(dx), abs(dy));
                if (dx < 0) sqXmin = xc - sqSide;
                if (dy < 0) sqYmin = yc - sqSide;
                HPEN hPenS = CreatePen(PS_SOLID, 1, currentColor);
                HGDIOBJ oldPenS = SelectObject(hdc, hPenS);
                Rectangle(hdc, sqXmin, sqYmin, sqXmin + sqSide, sqYmin + sqSide);
                SelectObject(hdc, oldPenS);
                DeleteObject(hPenS);
                shapes.push_back({ "ClipSquareWindow", sqXmin, sqYmin, sqSide, 0, 0, 0, 0, currentColor });
                cout << "Square clipping window: top-left(" << sqXmin << "," << sqYmin << "), side=" << sqSide << "\n";
                break;
            }
            case ID_CLIP_SQUARE_POINT: {
                bool inside = (xc >= sqXmin && xc <= sqXmin + sqSide && yc >= sqYmin && yc <= sqYmin + sqSide);
                SetPixel(hdc, xc, yc, inside ? RGB(0,0,255) : RGB(255,0,0));
                shapes.push_back({ "ClipSquarePoint", xc, yc, 0, 0, 0, 0, 0, inside ? RGB(0,0,255) : RGB(255,0,0) });
                cout << "Square clip point (" << xc << "," << yc << ") is " << (inside ? "inside\n" : "outside\n");
                break;
            }
            case ID_CLIP_SQUARE_LINE: {
                int cx1, cy1, cx2, cy2;
                int xmin = sqXmin, ymin = sqYmin;
                int xmax = sqXmin + sqSide, ymax = sqYmin + sqSide;
                if (cohenSutherlandClip(xc, yc, xe, ye, xmin, ymin, xmax, ymax, cx1, cy1, cx2, cy2)) {
                    lineBresenham(hdc, cx1, cy1, cx2, cy2, currentColor);
                    shapes.push_back({ "ClipSquareLine", cx1, cy1, cx2, cy2, sqXmin, sqYmin, sqSide, currentColor });
                    cout << "Clipped line in square window: (" << cx1 << "," << cy1 << ") -> (" << cx2 << "," << cy2 << ")\n";
                } else {
                    cout << "Line entirely outside square window\n";
                }
                break;
            }

            // Circular clipping (bonus)
            case ID_CLIP_CIRCLE_WINDOW: {
                windowX = xc;
                windowY = yc;
                R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
                windowR = R;
                circlePolar(hdc, xc, yc, R, currentColor);
                shapes.push_back({ "ClipCircleWindow", xc, yc, 0, 0, R, 0, 0, currentColor });
                cout << "Defined circular clipping window: center(" << xc << "," << yc << "), r=" << R << "\n";
                break;
            }
            case ID_CLIP_CIRCLE_LINE: {
                circleClipLine(hdc, xc, yc, xe, ye, windowX, windowY, windowR, currentColor);
                break;
            }
            case ID_CLIP_CIRCLE_POINT: {
                pointClipCircle(hdc, xc, yc, windowR, windowX, windowY, currentColor);
                shapes.push_back({ "ClipCirclePoint", xc, yc, 0, 0, windowX, windowY, windowR, currentColor });
                cout << "Circular clip test for point (" << xc << "," << yc << ")\n";
                break;
            }

            default:
                break;
        }
        break;

    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdcPaint = BeginPaint(hwnd, &ps);

        for (auto &s : shapes) {
            if (s.type == "DDA") {
                DDALine(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "MidpointLine") {
                lineBresenham(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "ParametricLine") {
                parametricLine(hdcPaint, s.x1, s.y1, s.x2, s.y2, 0, 0, 0, s.color);
            }
            else if (s.type == "CircleDirect") {
                circleDirectMethod(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CirclePolar") {
                circlePolar(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CircleIterative") {
                circleIterative(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CircleMidpoint") {
                midpointCircle(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "CircleFast") {
                circleFastBresenham(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "FillCircleLines") {
                fillWithLines(hdcPaint, s.x1, s.y1, 0, 0, s.a, s.color, s.b);
            }
            else if (s.type == "FillCircleCircles") {
                fillWithCircles(hdcPaint, s.x1, s.y1, 0, 0, s.a, s.color, s.b);
            }
            else if (s.type == "HermiteSquare") {
                POINT tl = { s.x1, s.y1 };
                FillHermiteSquare(hdcPaint, tl, s.a, s.color);
            }
            else if (s.type == "FillBezierRect") {
                POINT tl = { s.x1, s.y1 };
                POINT br = { s.x2, s.y2 };
                FillBezierRectangle(hdcPaint, tl, br, s.color);
            }
            else if (s.type == "ConvexPoly") {
                fillConvexPolygon(hdcPaint, s.verts, s.color);
            }
            else if (s.type == "NonConvexPoly") {
                fillNonConvexPolygon(hdcPaint, s.verts, s.color);
            }
            else if (s.type == "CardinalSpline") {
                drawCardinalSpline(hdcPaint, s.verts, s.color);
            }
            else if (s.type == "FloodRec") {
                COLORREF tar = GetPixel(hdcPaint, s.x1, s.y1);
                floodFillRec(hdcPaint, s.x1, s.y1, tar, s.color);
            }
            else if (s.type == "FloodStack") {
                COLORREF tar = GetPixel(hdcPaint, s.x1, s.y1);
                floodFillStack(hdcPaint, s.x1, s.y1, tar, s.color);
            }
            else if (s.type == "EllipseDirect") {
                ellipseDirect(hdcPaint, s.x1, s.y1, s.a, s.b, s.color);
            }
            else if (s.type == "EllipsePolar") {
                ellipsePolar(hdcPaint, s.x1, s.y1, s.a, s.b, s.color);
            }
            else if (s.type == "EllipseMidpoint") {
                ellipseMidpoint(hdcPaint, s.x1, s.y1, s.a, s.b, s.color);
            }
            else if (s.type == "ClipRectWindow") {
                drawRectangleWindow(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "ClipRectPoint") {
                SetPixel(hdcPaint, s.x1, s.y1, s.color);
            }
            else if (s.type == "ClipRectLine") {
                lineBresenham(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "ClipRectPoly") {
                if (!s.verts.empty()) {
                    HPEN hPenC = CreatePen(PS_SOLID, 1, s.color);
                    HGDIOBJ oldPenC = SelectObject(hdcPaint, hPenC);
                    MoveToEx(hdcPaint, s.verts[0].x, s.verts[0].y, NULL);
                    for (size_t i = 1; i < s.verts.size(); i++) {
                        LineTo(hdcPaint, s.verts[i].x, s.verts[i].y);
                    }
                    LineTo(hdcPaint, s.verts[0].x, s.verts[0].y);
                    SelectObject(hdcPaint, oldPenC);
                    DeleteObject(hPenC);
                }
            }
            else if (s.type == "ClipSquareWindow") {
                HPEN hPenS = CreatePen(PS_SOLID, 1, s.color);
                HGDIOBJ oldPenS = SelectObject(hdcPaint, hPenS);
                Rectangle(hdcPaint, s.x1, s.y1, s.x1 + s.x2, s.y1 + s.x2);
                SelectObject(hdcPaint, oldPenS);
                DeleteObject(hPenS);
            }
            else if (s.type == "ClipSquarePoint") {
                SetPixel(hdcPaint, s.x1, s.y1, s.color);
            }
            else if (s.type == "ClipSquareLine") {
                lineBresenham(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "ClipCircleWindow") {
                circlePolar(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "ClipCirclePoint") {
                SetPixel(hdcPaint, s.x1, s.y1, s.color);
            }
            else if (s.type == "ClipCircleLine") {
                lineBresenham(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
        }

        EndPaint(hwnd, &ps);
        break;
    }

    case WM_CREATE:
        CreateMenus(hwnd);
        break;

    case WM_DESTROY:
        PostQuitMessage(0);
        break;

    default:
        return DefWindowProc(hwnd, message, wParam, lParam);
    }

    return 0;
}

//-----------------------------------------------------------------------------
// WinMain: Entry Point
//-----------------------------------------------------------------------------

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    // Combine console and window
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    ShowWindow(GetConsoleWindow(), SW_SHOW);

    WNDCLASSEX wc = {};
    wc.cbSize        = sizeof(WNDCLASSEX);
    wc.style         = CS_DBLCLKS;
    wc.lpfnWndProc   = WindowProcedure;
    wc.hInstance     = hInstance;
    wc.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wc.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor       = LoadCursor(NULL, IDC_CROSS);          // crosshair cursor
    wc.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);  // white background
    wc.lpszClassName = _T("Clean2DDrawingApp");
    wc.lpszMenuName  = NULL;

    if (!RegisterClassEx(&wc)) return 0;

    HWND hwnd = CreateWindowEx(
        0,
        _T("Clean2DDrawingApp"),
        _T("2D Drawing Program"),
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 800, 600,
        HWND_DESKTOP, NULL, hInstance, NULL
    );

    ShowWindow(hwnd, nCmdShow);

    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return (int)msg.wParam;
}

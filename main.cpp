#include <tchar.h>
#include <windows.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stack>
#include <algorithm>

using namespace std;

// Forward declare the menu creation function
void menu(HWND hwnd);

//-------------------------------------------------------------------------------------------------
// Data Structures and Global Variables
//-------------------------------------------------------------------------------------------------

struct Shape {
    std::string type;           // e.g. "DDALine", "CirclePolar", "ClipRectLine", "ClipSquareLine", etc.
    int x1, y1;                 // Primary coordinates
    int x2, y2;                 // Secondary coordinates
    int a, b;                   // Extra parameters (radius, side length, etc.)
    int c;                      // Rarely used
    COLORREF color;             // Drawing color
    std::vector<POINT> vertices; // For polygons or splines
};

static std::vector<Shape> shapes;

static int xc, yc;      // First click
static int xe, ye;      // Second click
static int xe2, ye2;    // Third click (for polar ellipse)
static int R, R2;       // Radii or aux

static int windowX, windowY, windowR;     // Circular window
static int rectXmin, rectYmin, rectXmax, rectYmax; // Rectangular window
static int sqXmin, sqYmin, sqSide;        // Square clipping window

static std::vector<POINT> currentPolygon;
static std::vector<POINT> currentSplinePoints;

static HMENU hMenu;

// Mode IDs
static const int ID_FLOOD_REC            = 70;
static const int ID_FLOOD_STACK          = 71;
static const int ID_FILL_BEZIERRECT      = 80;
static const int ID_FILL_CONVEX_POLY     = 90;
static const int ID_FILL_NONCONVEX_POLY  = 91;
static const int ID_CARDINAL_SPLINE      = 95;
static const int ID_CLIP_RECT_WINDOW     = 100;
static const int ID_CLIP_RECT_POINT      = 101;
static const int ID_CLIP_RECT_LINE       = 102;
static const int ID_CLIP_RECT_POLY       = 103;
static const int ID_CLIP_SQUARE_WINDOW   = 110;
static const int ID_CLIP_SQUARE_POINT    = 111;
static const int ID_CLIP_SQUARE_LINE     = 112;

//-------------------------------------------------------------------------------------------------
// Utility Functions (LERP, Bezier, etc.)
//-------------------------------------------------------------------------------------------------

double LERP(double start, double end, double t) {
    return ((1 - t) * start) + (t * end);
}

double quadraticBezier(double x, double y, double z, double t) {
    double a = LERP(x, y, t);
    double b = LERP(y, z, t);
    return LERP(a, b, t);
}

double cubicBezier(double w, double x, double y, double z, double t) {
    double a = LERP(w, x, t);
    double b = LERP(x, y, t);
    double c = LERP(y, z, t);
    double d = LERP(a, b, t);
    double e = LERP(b, c, t);
    return LERP(d, e, t);
}

//-------------------------------------------------------------------------------------------------
// Hermite Curve Functions
//-------------------------------------------------------------------------------------------------

POINT EvaluateHermite(float t, POINT p0, POINT p1, POINT r0, POINT r1) {
    float h1 =  2 * t * t * t - 3 * t * t + 1;
    float h2 = -2 * t * t * t + 3 * t * t;
    float h3 =  t * t * t - 2 * t * t + t;
    float h4 =  t * t * t - t * t;
    POINT p;
    p.x = (int)round(h1 * p0.x + h2 * p1.x + h3 * r0.x + h4 * r1.x);
    p.y = (int)round(h1 * p0.y + h2 * p1.y + h3 * r0.y + h4 * r1.y);
    return p;
}

void FillHermiteSquare(HDC hdc, POINT topLeft, int size, COLORREF color) {
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    const int count = 50;
    for (int i = 0; i <= count; i++) {
        float t = (float)i / count;
        int x = topLeft.x + (int)round(t * size);
        POINT p0 = { x, topLeft.y };
        POINT p1 = { x, topLeft.y + size };
        POINT r0 = { (int)round(50.0f * sinf(t * 3.1415f)), 100 };
        POINT r1 = { (int)round(-50.0f * sinf(t * 3.1415f)), -100 };

        const int steps = 100;
        POINT prev = EvaluateHermite(0.0f, p0, p1, r0, r1);
        for (int j = 1; j <= steps; j++) {
            float tj = (float)j / steps;
            POINT curr = EvaluateHermite(tj, p0, p1, r0, r1);
            MoveToEx(hdc, prev.x, prev.y, NULL);
            LineTo(hdc, curr.x, curr.y);
            prev = curr;
        }
    }

    MoveToEx(hdc, topLeft.x, topLeft.y, NULL);
    LineTo(hdc, topLeft.x + size, topLeft.y);
    LineTo(hdc, topLeft.x + size, topLeft.y + size);
    LineTo(hdc, topLeft.x, topLeft.y + size);
    LineTo(hdc, topLeft.x, topLeft.y);

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

void drawCardinalSpline(HDC hdc, const std::vector<POINT>& pts, COLORREF color) {
    size_t n = pts.size();
    if (n < 3) return;
    std::vector<POINT> tangents(n);
    tangents[0] = {0, 0};
    tangents[n-1] = {0, 0};
    for (size_t i = 1; i < n - 1; ++i) {
        tangents[i].x = (int)round(0.5 * (pts[i + 1].x - pts[i - 1].x));
        tangents[i].y = (int)round(0.5 * (pts[i + 1].y - pts[i - 1].y));
    }

    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    const int stepsPerSegment = 50;
    for (size_t i = 0; i < n - 1; ++i) {
        POINT P0 = pts[i];
        POINT P1 = pts[i + 1];
        POINT R0 = tangents[i];
        POINT R1 = tangents[i + 1];

        POINT prev = EvaluateHermite(0.0f, P0, P1, R0, R1);
        for (int s = 1; s <= stepsPerSegment; ++s) {
            float t = (float)s / stepsPerSegment;
            POINT curr = EvaluateHermite(t, P0, P1, R0, R1);
            MoveToEx(hdc, prev.x, prev.y, NULL);
            LineTo(hdc, curr.x, curr.y);
            prev = curr;
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

//-------------------------------------------------------------------------------------------------
// Flood Fill Algorithms
//-------------------------------------------------------------------------------------------------

void floodFillRec(HDC hdc, int x, int y, COLORREF targetColor, COLORREF fillColor) {
    COLORREF current = GetPixel(hdc, x, y);
    if (current != targetColor || current == fillColor) return;
    SetPixel(hdc, x, y, fillColor);
    floodFillRec(hdc, x + 1, y, targetColor, fillColor);
    floodFillRec(hdc, x - 1, y, targetColor, fillColor);
    floodFillRec(hdc, x, y + 1, targetColor, fillColor);
    floodFillRec(hdc, x, y - 1, targetColor, fillColor);
}

void floodFillStack(HDC hdc, int x, int y, COLORREF targetColor, COLORREF fillColor) {
    if (targetColor == fillColor) return;
    std::stack<POINT> stk;
    stk.push({x, y});
    while (!stk.empty()) {
        POINT p = stk.top(); stk.pop();
        COLORREF current = GetPixel(hdc, p.x, p.y);
        if (current != targetColor || current == fillColor) continue;
        SetPixel(hdc, p.x, p.y, fillColor);
        stk.push({p.x + 1, p.y});
        stk.push({p.x - 1, p.y});
        stk.push({p.x, p.y + 1});
        stk.push({p.x, p.y - 1});
    }
}

//-------------------------------------------------------------------------------------------------
// Fill Rectangle with Horizontal Bezier Curves
//-------------------------------------------------------------------------------------------------

void FillBezierRectangle(HDC hdc, POINT rectTopLeft, POINT rectBotRight, COLORREF color) {
    int x1 = rectTopLeft.x, y1 = rectTopLeft.y;
    int x2 = rectBotRight.x, y2 = rectBotRight.y;
    if (x2 < x1) swap(x1, x2);
    if (y2 < y1) swap(y1, y2);

    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    int cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    const int segments = 50;

    for (int y = y1; y <= y2; y++) {
        double px0 = x1, py0 = y;
        double px1 = cx, py1 = cy;
        double px2 = x2, py2 = y;
        double prevX = px0, prevY = py0;
        for (int i = 1; i <= segments; i++) {
            double t = (double)i / segments;
            double xt = quadraticBezier(px0, px1, px2, t);
            double yt = quadraticBezier(py0, py1, py2, t);
            MoveToEx(hdc, (int)round(prevX), (int)round(prevY), NULL);
            LineTo(hdc, (int)round(xt), (int)round(yt));
            prevX = xt; prevY = yt;
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

//-------------------------------------------------------------------------------------------------
// Polygon Fill Algorithms (Convex & Non Convex)
//-------------------------------------------------------------------------------------------------

void fillConvexPolygon(HDC hdc, const std::vector<POINT>& verts, COLORREF color) {
    if (verts.size() < 3) return;
    int ymin = verts[0].y, ymax = verts[0].y;
    for (auto &p : verts) {
        if (p.y < ymin) ymin = p.y;
        if (p.y > ymax) ymax = p.y;
    }
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    for (int y = ymin; y <= ymax; ++y) {
        double leftX = 1e9, rightX = -1e9;
        for (size_t i = 0; i < verts.size(); ++i) {
            POINT p1 = verts[i], p2 = verts[(i + 1) % verts.size()];
            if (p1.y == p2.y) continue;
            int edgeYmin = (p1.y < p2.y ? p1.y : p2.y);
            int edgeYmax = (p1.y > p2.y ? p1.y : p2.y);
            if (y < edgeYmin || y > edgeYmax) continue;
            double x = p1.x + (double)(y - p1.y) * (p2.x - p1.x) / (double)(p2.y - p1.y);
            if (x < leftX)  leftX  = x;
            if (x > rightX) rightX = x;
        }
        if (leftX <= rightX) {
            MoveToEx(hdc, (int)ceil(leftX), y, NULL);
            LineTo(hdc, (int)floor(rightX), y);
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

void fillNonConvexPolygon(HDC hdc, const std::vector<POINT>& verts, COLORREF color) {
    if (verts.size() < 3) return;
    int ymin = verts[0].y, ymax = verts[0].y;
    for (auto &p : verts) {
        if (p.y < ymin) ymin = p.y;
        if (p.y > ymax) ymax = p.y;
    }
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HGDIOBJ oldPen = SelectObject(hdc, hPen);

    for (int y = ymin; y <= ymax; ++y) {
        vector<double> xIntersections;
        for (size_t i = 0; i < verts.size(); ++i) {
            POINT p1 = verts[i], p2 = verts[(i + 1) % verts.size()];
            if (p1.y == p2.y) continue;
            int edgeYmin = (p1.y < p2.y ? p1.y : p2.y);
            int edgeYmax = (p1.y > p2.y ? p1.y : p2.y);
            if (y < edgeYmin || y >= edgeYmax) continue;
            double x = p1.x + (double)(y - p1.y) * (p2.x - p1.x) / (double)(p2.y - p1.y);
            xIntersections.push_back(x);
        }
        sort(xIntersections.begin(), xIntersections.end());
        for (size_t i = 0; i + 1 < xIntersections.size(); i += 2) {
            int xStart = (int)ceil(xIntersections[i]);
            int xEnd   = (int)floor(xIntersections[i + 1]);
            if (xStart <= xEnd) {
                MoveToEx(hdc, xStart, y, NULL);
                LineTo(hdc, xEnd,   y);
            }
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(hPen);
}

//-------------------------------------------------------------------------------------------------
// Rectangle & Square Clipping Utility Functions
//-------------------------------------------------------------------------------------------------

void drawRectangleWindow(HDC hdc, int xmin, int ymin, int xmax, int ymax, COLORREF color) {
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

// Cohen–Sutherland outcode computation
const int INSIDE_CS = 0; // 0000
const int LEFT_CS   = 1; // 0001
const int RIGHT_CS  = 2; // 0010
const int BOTTOM_CS = 4; // 0100
const int TOP_CS    = 8; // 1000

int computeOutCode(int x, int y, int xmin, int ymin, int xmax, int ymax) {
    int code = INSIDE_CS;
    if (x < xmin)       code |= LEFT_CS;
    else if (x > xmax)  code |= RIGHT_CS;
    if (y < ymin)       code |= BOTTOM_CS;
    else if (y > ymax)  code |= TOP_CS;
    return code;
}

bool cohenSutherlandClip(int x1, int y1, int x2, int y2,
                         int xmin, int ymin, int xmax, int ymax,
                         int &cx1, int &cy1, int &cx2, int &cy2)
{
    int outcode1 = computeOutCode(x1, y1, xmin, ymin, xmax, ymax);
    int outcode2 = computeOutCode(x2, y2, xmin, ymin, xmax, ymax);
    bool accept = false;
    double xd = x2 - x1, yd = y2 - y1;

    while (true) {
        if ((outcode1 | outcode2) == 0) {
            accept = true;
            cx1 = x1; cy1 = y1;
            cx2 = x2; cy2 = y2;
            break;
        }
        else if (outcode1 & outcode2) {
            break;
        }
        else {
            int outcodeOut = outcode1 ? outcode1 : outcode2;
            double x, y;
            if (outcodeOut & TOP_CS) {
                x = x1 + xd * (ymax - y1) / yd;  y = ymax;
            }
            else if (outcodeOut & BOTTOM_CS) {
                x = x1 + xd * (ymin - y1) / yd;  y = ymin;
            }
            else if (outcodeOut & RIGHT_CS) {
                y = y1 + yd * (xmax - x1) / xd;  x = xmax;
            }
            else {
                y = y1 + yd * (xmin - x1) / xd;  x = xmin;
            }

            if (outcodeOut == outcode1) {
                x1 = (int)round(x); y1 = (int)round(y);
                outcode1 = computeOutCode(x1, y1, xmin, ymin, xmax, ymax);
            }
            else {
                x2 = (int)round(x); y2 = (int)round(y);
                outcode2 = computeOutCode(x2, y2, xmin, ymin, xmax, ymax);
            }
        }
    }
    return accept;
}

// Sutherland–Hodgman polygon clipping
inline POINT intersectEdge(const POINT &P1, const POINT &P2, int edge,
                           int xmin, int ymin, int xmax, int ymax)
{
    double x1 = P1.x, y1 = P1.y;
    double x2 = P2.x, y2 = P2.y;
    double dx = x2 - x1, dy = y2 - y1;
    double x, y;

    switch (edge) {
        case 0: // LEFT: x = xmin
            x = xmin;
            y = y1 + dy * (xmin - x1) / dx;
            break;
        case 1: // RIGHT: x = xmax
            x = xmax;
            y = y1 + dy * (xmax - x1) / dx;
            break;
        case 2: // BOTTOM: y = ymin
            y = ymin;
            x = x1 + dx * (ymin - y1) / dy;
            break;
        case 3: // TOP: y = ymax
            y = ymax;
            x = x1 + dx * (ymax - y1) / dy;
            break;
        default:
            x = x1; y = y1;
    }
    return { (int)round(x), (int)round(y) };
}

inline bool insideEdge(const POINT &P, int edge, int xmin, int ymin, int xmax, int ymax) {
    switch (edge) {
        case 0: return P.x >= xmin;  // LEFT
        case 1: return P.x <= xmax;  // RIGHT
        case 2: return P.y >= ymin;  // BOTTOM
        case 3: return P.y <= ymax;  // TOP
    }
    return false;
}

void sutherlandHodgmanPolygonClip(const std::vector<POINT> &inPoly,
                                  std::vector<POINT> &outPoly,
                                  int xmin, int ymin, int xmax, int ymax)
{
    std::vector<POINT> tempPoly = inPoly;

    for (int edge = 0; edge < 4; ++edge) {
        outPoly.clear();
        if (tempPoly.empty()) break;
        POINT S = tempPoly.back();

        for (const POINT &E : tempPoly) {
            bool insideE = insideEdge(E, edge, xmin, ymin, xmax, ymax);
            bool insideS = insideEdge(S, edge, xmin, ymin, xmax, ymax);

            if (insideE) {
                if (!insideS) {
                    POINT I = intersectEdge(S, E, edge, xmin, ymin, xmax, ymax);
                    outPoly.push_back(I);
                }
                outPoly.push_back(E);
            }
            else if (insideS) {
                POINT I = intersectEdge(S, E, edge, xmin, ymin, xmax, ymax);
                outPoly.push_back(I);
            }
            S = E;
        }
        tempPoly = outPoly;
    }
}

//-------------------------------------------------------------------------------------------------
// Line Drawing Algorithms
//-------------------------------------------------------------------------------------------------

void lineBresenham(HDC hdc, int x0, int y0, int x1, int y1, COLORREF color) {
    int dx = abs(x1 - x0), dy = abs(y1 - y0);
    int x = x0, y = y0;
    int sx = (x0 < x1 ? 1 : -1);
    int sy = (y0 < y1 ? 1 : -1);
    int err = dx - dy;

    while (true) {
        SetPixel(hdc, x, y, color);
        if (x == x1 && y == y1) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x += sx; }
        if (e2 < dx)  { err += dx; y += sy; }
    }
}

void DDALine(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    int dx = x2 - x1, dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    float Xinc = dx / (float)steps;
    float Yinc = dy / (float)steps;
    float X = (float)x1, Y = (float)y1;
    for (int i = 0; i <= steps; i++) {
        SetPixel(hdc, (int)round(X), (int)round(Y), color);
        X += Xinc; Y += Yinc;
    }
}

void parametricLine(HDC hdc, int x1, int y1, int x2, int y2,
                    int cx, int cy, int radius, COLORREF color) {
    int dx = x2 - x1, dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    if (steps == 0) {
        if (!radius || ((x1 - cx)*(x1 - cx) + (y1 - cy)*(y1 - cy) <= radius*radius))
            SetPixel(hdc, x1, y1, color);
        return;
    }
    double tStep = 1.0 / steps;
    for (double t = 0.0; t <= 1.0; t += tStep) {
        int x = (int)round(LERP(x1, x2, t));
        int y = (int)round(LERP(y1, y2, t));
        if (!radius || ((x - cx)*(x - cx) + (y - cy)*(y - cy) <= radius*radius))
            SetPixel(hdc, x, y, color);
    }
}

//-------------------------------------------------------------------------------------------------
// Circle Drawing Algorithms
//-------------------------------------------------------------------------------------------------

void draw8Points(HDC hdc, int cx, int cy, int x, int y, COLORREF color) {
    SetPixel(hdc, cx + x, cy + y, color);
    SetPixel(hdc, cx - x, cy + y, color);
    SetPixel(hdc, cx - x, cy - y, color);
    SetPixel(hdc, cx + x, cy - y, color);
    SetPixel(hdc, cx + y, cy + x, color);
    SetPixel(hdc, cx - y, cy + x, color);
    SetPixel(hdc, cx - y, cy - x, color);
    SetPixel(hdc, cx + y, cy - x, color);
}

void midpointCircle(HDC hdc, int cx, int cy, int radius, COLORREF color) {
    int x = 0, y = radius, d = 1 - radius;
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

void circleDirectMethod(HDC hdc, int cx, int cy, int R, COLORREF color) {
    double x = 0.0, y = (double)R, R2 = (double)R * R;
    while (x < y) {
        draw8Points(hdc, cx, cy, (int)round(x), (int)round(y), color);
        x += 0.1;
        y = sqrt(R2 - x * x);
    }
}

void circlePolar(HDC hdc, int cx, int cy, int R, COLORREF color) {
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

void circleIterative(HDC hdc, int cx, int cy, int R, COLORREF color) {
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

void circleFastBresenham(HDC hdc, int cx, int cy, int R, COLORREF color) {
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

void drawQuadrant(HDC hdc, int cx, int cy, int x, int y, int quadrant, COLORREF color) {
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
    }
}

void circleFastBresenhamQuadrant(HDC hdc, int cx, int cy, int R, int quadrant, COLORREF color) {
    int x = 0, y = R;
    int d = 1 - R, c1 = 3, c2 = 5 - 2 * R;
    drawQuadrant(hdc, cx, cy, x, y, quadrant, color);
    while (x < y) {
        if (d < 0) {
            d += c1; c2 += 2;
        } else {
            d += c2; c2 += 4; y--;
        }
        x++; c1 += 2;
        drawQuadrant(hdc, cx, cy, x, y, quadrant, color);
    }
}

//-------------------------------------------------------------------------------------------------
// Ellipse Drawing Algorithms
//-------------------------------------------------------------------------------------------------

void draw4Points(HDC hdc, int cx, int cy, int x, int y, COLORREF color) {
    SetPixel(hdc, cx + x, cy + y, color);
    SetPixel(hdc, cx - x, cy + y, color);
    SetPixel(hdc, cx - x, cy - y, color);
    SetPixel(hdc, cx + x, cy - y, color);
}

void ellipseDirect(HDC hdc, int cx, int cy, int a, int b, COLORREF color) {
    int x = 0, y = b;
    draw4Points(hdc, cx, cy, x, y, color);
    while (x * (b*b) < (a*a) * y) {
        x++;
        y = (int)round(sqrt((double)(b*b*(a*a - x*x)) / (a*a)));
        draw4Points(hdc, cx, cy, x, y, color);
    }
    y = 0; x = a;
    draw4Points(hdc, cx, cy, x, y, color);
    while (x * (b*b) > (a*a) * y) {
        y++;
        x = (int)round(sqrt((double)(a*a*(b*b - y*y)) / (b*b)));
        draw4Points(hdc, cx, cy, x, y, color);
    }
}

void ellipsePolar(HDC hdc, int cx, int cy, int R1, int R2, COLORREF color) {
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

//-------------------------------------------------------------------------------------------------
// Filling Algorithms (Lines & Circles for quarter circle)
//-------------------------------------------------------------------------------------------------

void fillWithLines(HDC hdc, int cx, int cy,
                   int unusedA, int unusedB, int R, COLORREF color, int q)
{
    int x = R, y = 0;
    double theta = 0.0, dtheta = 1.0 / R;
    switch (q) {
        case 1:  lineBresenham(hdc, cx, cy, cx + x, cy - y, color); break;
        case 2:  lineBresenham(hdc, cx, cy, cx + x, cy + y, color); break;
        case 3:  lineBresenham(hdc, cx, cy, cx - x, cy + y, color); break;
        case 4:  lineBresenham(hdc, cx, cy, cx - x, cy - y, color); break;
    }
    while (x * R > y) {
        theta += dtheta;
        x = (int)round(R * cos(theta));
        y = (int)round(R * sin(theta));
        switch (q) {
            case 1: lineBresenham(hdc, cx, cy, cx + x, cy - y, color); break;
            case 2: lineBresenham(hdc, cx, cy, cx + x, cy + y, color); break;
            case 3: lineBresenham(hdc, cx, cy, cx - x, cy + y, color); break;
            case 4: lineBresenham(hdc, cx, cy, cx - x, cy - y, color); break;
        }
    }
}

void fillWithCircles(HDC hdc, int cx, int cy,
                     int unusedA, int unusedB, int R, COLORREF color, int q)
{
    int circleIndex = 0;
    switch (q) {
        case 1: circleIndex = 0; break;
        case 2: circleIndex = 3; break;
        case 3: circleIndex = 2; break;
        case 4: circleIndex = 1; break;
    }
    for (int r = 0; r < R; r++) {
        circleFastBresenhamQuadrant(hdc, cx, cy, r, circleIndex, color);
    }
}

//-------------------------------------------------------------------------------------------------
// Point Clipping inside a Circle
//-------------------------------------------------------------------------------------------------

void pointClip(HDC hdc, int x, int y,
               int radius, int centerX, int centerY, COLORREF color)
{
    double dx2 = pow(centerX - x, 2);
    double dy2 = pow(centerY - y, 2);
    double dist = sqrt(dx2 + dy2);
    if (dist <= radius)
        SetPixel(hdc, x, y, RGB(0, 0, 255));
    else
        SetPixel(hdc, x, y, RGB(255, 0, 0));
}

//-------------------------------------------------------------------------------------------------
// File I/O: Save & Load all shapes
//-------------------------------------------------------------------------------------------------

void PerformSave(const char* filename) {
    ofstream out(filename);
    if (!out.is_open()) return;

    for (const auto& s : shapes) {
        out << s.type << " "
            << s.x1 << " " << s.y1 << " "
            << s.x2 << " " << s.y2 << " "
            << s.a  << " " << s.b  << " "
            << s.c  << " "
            << s.color;
        if (!s.vertices.empty()) {
            out << " " << s.vertices.size();
            for (auto &p : s.vertices) {
                out << " " << p.x << " " << p.y;
            }
        }
        out << "\n";
    }
    out.close();
}

void PerformLoad(const char* filename) {
    ifstream in(filename);
    if (!in.is_open()) return;

    shapes.clear();
    string line;
    while (std::getline(in, line)) {
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
            s.type == "ClipSquarePoly")
        {
            size_t n;
            iss >> n;
            s.vertices.resize(n);
            for (size_t i = 0; i < n; i++) {
                iss >> s.vertices[i].x >> s.vertices[i].y;
            }
        }
        shapes.push_back(s);
    }
    in.close();
    InvalidateRect(GetActiveWindow(), NULL, TRUE);
}

//-------------------------------------------------------------------------------------------------
// Window Procedure & Message Handling
//-------------------------------------------------------------------------------------------------

LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    static int currentShape = 0;
    static COLORREF currentColor = RGB(0, 0, 0);

    HDC hdc = GetDC(hwnd);

    switch (message) {
    case WM_COMMAND:
        if (wParam == 1) {
            PerformSave("shapes_data.txt");
            ReleaseDC(hwnd, hdc);
            cout << "Action: Save shapes\n";
        }
        else if (wParam == 2) {
            PerformLoad("shapes_data.txt");
            cout << "Action: Load shapes\n";
        }
        else if (wParam == 3) {
            shapes.clear();
            currentPolygon.clear();
            currentSplinePoints.clear();
            InvalidateRect(hwnd, NULL, TRUE);
            cout << "Action: Clear canvas\n";
        }
        else if (wParam == 41) { currentColor = RGB(255, 0, 0); cout << "Color: Red\n"; }
        else if (wParam == 42) { currentColor = RGB(0, 0, 0); cout << "Color: Black\n"; }
        else if (wParam == 43) { currentColor = RGB(0, 0, 255); cout << "Color: Blue\n"; }
        else if (wParam == 44) { currentColor = RGB(255, 255, 255); cout << "Color: White\n"; }
        else if (wParam == 45) { currentColor = RGB(255, 63, 127); cout << "Color: Pink\n"; }

        else if (wParam == ID_FLOOD_REC) {
            currentShape = ID_FLOOD_REC;
            cout << "Mode: Flood Fill (Recursive)\n";
        }
        else if (wParam == ID_FLOOD_STACK) {
            currentShape = ID_FLOOD_STACK;
            cout << "Mode: Flood Fill (Non-Recursive)\n";
        }
        else if (wParam == ID_FILL_BEZIERRECT) {
            currentShape = ID_FILL_BEZIERRECT;
            cout << "Mode: Fill Rectangle (Bezier)\n";
        }
        else if (wParam == ID_FILL_CONVEX_POLY) {
            currentShape = ID_FILL_CONVEX_POLY;
            currentPolygon.clear();
            cout << "Mode: Fill Convex Polygon (click vertices, right click to finish)\n";
        }
        else if (wParam == ID_FILL_NONCONVEX_POLY) {
            currentShape = ID_FILL_NONCONVEX_POLY;
            currentPolygon.clear();
            cout << "Mode: Fill Non Convex Polygon (click vertices, right click to finish)\n";
        }
        else if (wParam == ID_CARDINAL_SPLINE) {
            currentShape = ID_CARDINAL_SPLINE;
            currentSplinePoints.clear();
            cout << "Mode: Cardinal Spline (click control points, right click to finish)\n";
        }

        else if (wParam == 51) {
            currentShape = 51;
            cout << "Mode: Define Circular Clipping Window (click center)\n";
        }
        else if (wParam == 52) {
            currentShape = 52;
            cout << "Mode: Clip Line in Circular Window (click start)\n";
        }
        else if (wParam == 53) {
            currentShape = 53;
            cout << "Mode: Clip Point in Circular Window (click point, right click to clip)\n";
        }

        else if (wParam == ID_CLIP_RECT_WINDOW) {
            currentShape = ID_CLIP_RECT_WINDOW;
            cout << "Mode: Define Rectangular Clipping Window (click first corner)\n";
        }
        else if (wParam == ID_CLIP_RECT_POINT) {
            currentShape = ID_CLIP_RECT_POINT;
            cout << "Mode: Clip Point in Rectangular Window (click point, right click to clip)\n";
        }
        else if (wParam == ID_CLIP_RECT_LINE) {
            currentShape = ID_CLIP_RECT_LINE;
            cout << "Mode: Clip Line in Rectangular Window (click first corner of line)\n";
        }
        else if (wParam == ID_CLIP_RECT_POLY) {
            currentShape = ID_CLIP_RECT_POLY;
            currentPolygon.clear();
            cout << "Mode: Clip Polygon in Rectangular Window (click vertices, right click to finish)\n";
        }

        else if (wParam == ID_CLIP_SQUARE_WINDOW) {
            currentShape = ID_CLIP_SQUARE_WINDOW;
            cout << "Mode: Define Square Clipping Window (click first corner)\n";
        }
        else if (wParam == ID_CLIP_SQUARE_POINT) {
            currentShape = ID_CLIP_SQUARE_POINT;
            cout << "Mode: Clip Point in Square Window (click point, right click to clip)\n";
        }
        else if (wParam == ID_CLIP_SQUARE_LINE) {
            currentShape = ID_CLIP_SQUARE_LINE;
            cout << "Mode: Clip Line in Square Window (click first corner of line)\n";
        }

        else if (wParam >= 101 && wParam <= 104 && wParam < 200) {
            // Line-filling quadrants handled in WM_RBUTTONDOWN
        }
        else if (wParam >= 201 && wParam <= 204) {
            int quadrant = (int)wParam - 200;
            bool foundCircle = false;
            int cx = 0, cy = 0, radius = 0;
            for (int i = (int)shapes.size() - 1; i >= 0; --i) {
                const auto &s = shapes[i];
                if (s.type == "CircleDirect" ||
                    s.type == "CirclePolar"  ||
                    s.type == "CircleIter"   ||
                    s.type == "CircleMid"    ||
                    s.type == "CircleFast")
                {
                    foundCircle = true;
                    cx = s.x1; cy = s.y1; radius = s.a;
                    break;
                }
            }
            if (foundCircle) {
                fillWithCircles(hdc, cx, cy, 0, 0, radius, currentColor, quadrant);
                shapes.push_back({ "FillCircles", cx, cy, 0, 0, radius, quadrant, 0, currentColor });
                cout << "Filled last circle with circles (quadrant " << quadrant << ")\n";
            } else {
                cout << "No circle found to fill!\n";
            }
        }

        else {
            currentShape = (int)wParam;
            switch (currentShape) {
                case 11:  cout << "Mode: DDA Line\n"; break;
                case 12:  cout << "Mode: Midpoint Line\n"; break;
                case 13:  cout << "Mode: Parametric Line\n"; break;
                case 14:  cout << "Mode: Circle Direct\n"; break;
                case 15:  cout << "Mode: Circle Polar\n"; break;
                case 16:  cout << "Mode: Circle Iterative\n"; break;
                case 17:  cout << "Mode: Circle Midpoint\n"; break;
                case 18:  cout << "Mode: Circle Fast Bresenham\n"; break;
                case 19:  cout << "Mode: Ellipse Direct\n"; break;
                case 20:  cout << "Mode: Ellipse Polar\n"; break;
                case 51:  cout << "Mode: Define Circular Window\n"; break;
                case 52:  cout << "Mode: Clip Line in Circle\n"; break;
                case 53:  cout << "Mode: Clip Point in Circle\n"; break;
                case 60:  cout << "Mode: Hermite Square\n"; break;
                case ID_CLIP_SQUARE_WINDOW: cout << "Mode: Define Square Window\n"; break;
                case ID_CLIP_SQUARE_POINT:  cout << "Mode: Clip Point in Square\n"; break;
                case ID_CLIP_SQUARE_LINE:   cout << "Mode: Clip Line in Square\n"; break;
                default:  break;
            }
        }
        break;

    case WM_LBUTTONDOWN:
        xc = LOWORD(lParam);
        yc = HIWORD(lParam);

        if (currentShape == 51) {
            windowX = xc; windowY = yc;
            cout << "Set circular clip center: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == 53) {
            cout << "Selected circular clip point: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == 60) {
            cout << "Set Hermite square top-left: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_FLOOD_REC || currentShape == ID_FLOOD_STACK) {
            cout << "Flood fill seed: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_FILL_BEZIERRECT) {
            cout << "Set Bezier rectangle top-left: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_FILL_CONVEX_POLY || currentShape == ID_FILL_NONCONVEX_POLY) {
            currentPolygon.push_back({ xc, yc });
            SetPixel(hdc, xc, yc, currentColor);
            cout << "Added polygon vertex: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CARDINAL_SPLINE) {
            currentSplinePoints.push_back({ xc, yc });
            SetPixel(hdc, xc, yc, currentColor);
            cout << "Added spline control point: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CLIP_RECT_WINDOW) {
            rectXmin = xc; rectYmin = yc;
            cout << "Set rectangular clip window first corner: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CLIP_RECT_POINT) {
            cout << "Selected rectangular clip point: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CLIP_RECT_LINE) {
            cout << "Set rectangular clip line first point: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CLIP_RECT_POLY) {
            currentPolygon.push_back({ xc, yc });
            SetPixel(hdc, xc, yc, currentColor);
            cout << "Added clip polygon vertex: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CLIP_SQUARE_WINDOW) {
            sqXmin = xc; sqYmin = yc;
            cout << "Set square clip window first corner: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CLIP_SQUARE_POINT) {
            cout << "Selected square clip point: (" << xc << ", " << yc << ")\n";
        }
        else if (currentShape == ID_CLIP_SQUARE_LINE) {
            cout << "Set square clip line first point: (" << xc << ", " << yc << ")\n";
        }
        break;

    case WM_RBUTTONDOWN:
        xe = LOWORD(lParam);
        ye = HIWORD(lParam);

        if (currentShape == 11) {
            DDALine(hdc, xc, yc, xe, ye, currentColor);
            shapes.push_back({ "DDALine", xc, yc, xe, ye, 0, 0, 0, currentColor });
            cout << "Drew DDA line: (" << xc << "," << yc << ") -> (" << xe << "," << ye << ")\n";
        }
        else if (currentShape == 12) {
            lineBresenham(hdc, xc, yc, xe, ye, currentColor);
            shapes.push_back({ "MidLine", xc, yc, xe, ye, 0, 0, 0, currentColor });
            cout << "Drew Midpoint line: (" << xc << "," << yc << ") -> (" << xe << "," << ye << ")\n";
        }
        else if (currentShape == 13) {
            parametricLine(hdc, xc, yc, xe, ye, 0, 0, 0, currentColor);
            shapes.push_back({ "ParamLine", xc, yc, xe, ye, 0, 0, 0, currentColor });
            cout << "Drew Parametric line: (" << xc << "," << yc << ") -> (" << xe << "," << ye << ")\n";
        }
        else if (currentShape == 14) {
            R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circleDirectMethod(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleDirect", xc, yc, 0, 0, R, 0, 0, currentColor });
            cout << "Drew Direct circle: center (" << xc << "," << yc << "), r=" << R << "\n";
        }
        else if (currentShape == 15) {
            R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circlePolar(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CirclePolar", xc, yc, 0, 0, R, 0, 0, currentColor });
            cout << "Drew Polar circle: center (" << xc << "," << yc << "), r=" << R << "\n";
        }
        else if (currentShape == 16) {
            R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circleIterative(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleIter", xc, yc, 0, 0, R, 0, 0, currentColor });
            cout << "Drew Iterative circle: center (" << xc << "," << yc << "), r=" << R << "\n";
        }
        else if (currentShape == 17) {
            R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            midpointCircle(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleMid", xc, yc, 0, 0, R, 0, 0, currentColor });
            cout << "Drew Midpoint circle: center (" << xc << "," << yc << "), r=" << R << "\n";
        }
        else if (currentShape == 18) {
            R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            circleFastBresenham(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "CircleFast", xc, yc, 0, 0, R, 0, 0, currentColor });
            cout << "Drew FastBresenham circle: center (" << xc << "," << yc << "), r=" << R << "\n";
        }
        else if (currentShape == 19) {
            int a = abs(xe - xc), b = abs(ye - yc);
            ellipseDirect(hdc, xc, yc, a, b, currentColor);
            shapes.push_back({ "EllipseDirect", xc, yc, 0, 0, a, b, 0, currentColor });
            cout << "Drew Direct ellipse: center (" << xc << "," << yc << "), a=" << a << ", b=" << b << "\n";
        }
        else if (currentShape == 51) {
            windowX = xc; windowY = yc;
            R = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            windowR = R;
            circlePolar(hdc, xc, yc, R, currentColor);
            shapes.push_back({ "ClipWindow", xc, yc, 0, 0, R, 0, 0, currentColor });
            cout << "Defined circular clip window: center (" << xc << "," << yc << "), r=" << R << "\n";
        }
        else if (currentShape == 52) {
            int cx1, cy1, cx2, cy2;
            if (cohenSutherlandClip(xc, yc, xe, ye,
                                    windowX - windowR, windowY - windowR,
                                    windowX + windowR, windowY + windowR,
                                    cx1, cy1, cx2, cy2))
            {
                lineBresenham(hdc, cx1, cy1, cx2, cy2, currentColor);
                shapes.push_back({ "ClipLine", cx1, cy1, cx2, cy2, windowX, windowY, windowR, currentColor });
                cout << "Clipped line in circular window.\n";
            } else {
                cout << "Line entirely outside the circular window.\n";
            }
        }
        else if (currentShape == 53) {
            pointClip(hdc, xc, yc, windowR, windowX, windowY, currentColor);
            shapes.push_back({ "ClipPoint", xc, yc, 0, 0, windowX, windowY, windowR, currentColor });
            cout << "Did circular clipping test on point (" << xc << "," << yc << ")\n";
        }
        else if (currentShape == 60) {
            POINT topLeft = { xc, yc };
            int dx = abs(xe - xc), dy = abs(ye - yc);
            int size = min(dx, dy);
            FillHermiteSquare(hdc, topLeft, size, currentColor);
            shapes.push_back({ "HermiteSquare", topLeft.x, topLeft.y, 0, 0, size, 0, 0, currentColor });
            cout << "Drew Hermite square top-left (" << xc << "," << yc << "), size=" << size << "\n";
        }
        else if (currentShape == ID_FLOOD_REC || currentShape == ID_FLOOD_STACK) {
            COLORREF targetColor = GetPixel(hdc, xc, yc);
            if (targetColor != currentColor) {
                if (currentShape == ID_FLOOD_REC) {
                    floodFillRec(hdc, xc, yc, targetColor, currentColor);
                    shapes.push_back({ "FloodRec", xc, yc, 0, 0, 0, 0, 0, currentColor });
                    cout << "Performed recursive flood-fill at (" << xc << "," << yc << ")\n";
                } else {
                    floodFillStack(hdc, xc, yc, targetColor, currentColor);
                    shapes.push_back({ "FloodStack", xc, yc, 0, 0, 0, 0, 0, currentColor });
                    cout << "Performed stack-based flood-fill at (" << xc << "," << yc << ")\n";
                }
            } else {
                cout << "Seed color = fill color; nothing filled.\n";
            }
        }
        else if (currentShape == ID_FILL_BEZIERRECT) {
            int left   = min(xc, xe), right  = max(xc, xe);
            int top    = min(yc, ye), bottom = max(yc, ye);
            POINT topLeft     = { left, top };
            POINT bottomRight = { right, bottom };
            FillBezierRectangle(hdc, topLeft, bottomRight, currentColor);
            shapes.push_back({ "FillBezierRect", left, top, right, bottom, 0, 0, 0, currentColor });
            cout << "Filled rectangle with Bezier: top-left(" << left << "," << top << "), bottom-right("
                 << right << "," << bottom << ")\n";
        }
        else if (currentShape == ID_FILL_CONVEX_POLY) {
            if (currentPolygon.size() >= 3) {
                fillConvexPolygon(hdc, currentPolygon, currentColor);
                Shape polyShape;
                polyShape.type = "ConvexPoly";
                polyShape.x1 = polyShape.y1 = polyShape.x2 = polyShape.y2 = polyShape.a = polyShape.b = polyShape.c = 0;
                polyShape.color = currentColor;
                polyShape.vertices = currentPolygon;
                shapes.push_back(polyShape);
                cout << "Filled convex polygon with " << currentPolygon.size() << " vertices.\n";
            } else {
                cout << "Not enough vertices for convex polygon.\n";
            }
            currentPolygon.clear();
        }
        else if (currentShape == ID_FILL_NONCONVEX_POLY) {
            if (currentPolygon.size() >= 3) {
                fillNonConvexPolygon(hdc, currentPolygon, currentColor);
                Shape polyShape;
                polyShape.type = "NonConvexPoly";
                polyShape.x1 = polyShape.y1 = polyShape.x2 = polyShape.y2 = polyShape.a = polyShape.b = polyShape.c = 0;
                polyShape.color = currentColor;
                polyShape.vertices = currentPolygon;
                shapes.push_back(polyShape);
                cout << "Filled non-convex polygon with " << currentPolygon.size() << " vertices.\n";
            } else {
                cout << "Not enough vertices for non-convex polygon.\n";
            }
            currentPolygon.clear();
        }
        else if (currentShape == ID_CARDINAL_SPLINE) {
            if (currentSplinePoints.size() >= 3) {
                drawCardinalSpline(hdc, currentSplinePoints, currentColor);
                Shape splineShape;
                splineShape.type = "CardinalSpline";
                splineShape.x1 = splineShape.y1 = splineShape.x2 = splineShape.y2 = splineShape.a = splineShape.b = splineShape.c = 0;
                splineShape.color = currentColor;
                splineShape.vertices = currentSplinePoints;
                shapes.push_back(splineShape);
                cout << "Drew Cardinal spline with " << currentSplinePoints.size() << " control points.\n";
            } else {
                cout << "Not enough control points for spline.\n";
            }
            currentSplinePoints.clear();
        }
        else if (currentShape == ID_CLIP_RECT_WINDOW) {
            rectXmin = xc;  rectYmin = yc;
            rectXmax = xe;  rectYmax = ye;
            if (rectXmax < rectXmin) swap(rectXmin, rectXmax);
            if (rectYmax < rectYmin) swap(rectYmin, rectYmax);
            drawRectangleWindow(hdc, rectXmin, rectYmin, rectXmax, rectYmax, currentColor);
            shapes.push_back({ "ClipRectWindow", rectXmin, rectYmin, rectXmax, rectYmax, 0, 0, 0, currentColor });
            cout << "Defined rectangularClip window: (" << rectXmin << "," << rectYmin << ") to ("
                 << rectXmax << "," << rectYmax << ")\n";
        }
        else if (currentShape == ID_CLIP_RECT_POINT) {
            bool inside = (xc >= rectXmin && xc <= rectXmax && yc >= rectYmin && yc <= rectYmax);
            SetPixel(hdc, xc, yc, inside ? RGB(0, 0, 255) : RGB(255, 0, 0));
            shapes.push_back({ "ClipRectPoint", xc, yc, 0, 0, 0, 0, 0, inside ? RGB(0, 0, 255) : RGB(255, 0, 0) });
            cout << "Rectangular clip point (" << xc << "," << yc << ") is " << (inside ? "inside\n" : "outside\n");
        }
        else if (currentShape == ID_CLIP_RECT_LINE) {
            int cx1, cy1, cx2, cy2;
            if (cohenSutherlandClip(xc, yc, xe, ye, rectXmin, rectYmin, rectXmax, rectYmax, cx1, cy1, cx2, cy2)) {
                lineBresenham(hdc, cx1, cy1, cx2, cy2, currentColor);
                Shape clipLineShape;
                clipLineShape.type = "ClipRectLine";
                clipLineShape.x1 = cx1; clipLineShape.y1 = cy1;
                clipLineShape.x2 = cx2; clipLineShape.y2 = cy2;
                clipLineShape.a = rectXmin; clipLineShape.b = rectYmin;
                clipLineShape.c = 0; clipLineShape.color = currentColor;
                shapes.push_back(clipLineShape);
                cout << "Clipped line in rectangular window.\n";
            } else {
                cout << "Line entirely outside rectangular window.\n";
            }
        }
        else if (currentShape == ID_CLIP_RECT_POLY) {
            if (currentPolygon.size() >= 3) {
                vector<POINT> clippedPoly;
                sutherlandHodgmanPolygonClip(currentPolygon, clippedPoly,
                                             rectXmin, rectYmin, rectXmax, rectYmax);
                if (!clippedPoly.empty()) {
                    HPEN hPenC = CreatePen(PS_SOLID, 1, currentColor);
                    HGDIOBJ oldPenC = SelectObject(hdc, hPenC);
                    MoveToEx(hdc, clippedPoly[0].x, clippedPoly[0].y, NULL);
                    for (size_t i = 1; i < clippedPoly.size(); ++i) {
                        LineTo(hdc, clippedPoly[i].x, clippedPoly[i].y);
                    }
                    LineTo(hdc, clippedPoly[0].x, clippedPoly[0].y);
                    SelectObject(hdc, oldPenC);
                    DeleteObject(hPenC);

                    Shape clipPolyShape;
                    clipPolyShape.type = "ClipRectPoly";
                    clipPolyShape.x1 = clipPolyShape.y1 = clipPolyShape.x2 = clipPolyShape.y2 = 0;
                    clipPolyShape.a = clipPolyShape.b = clipPolyShape.c = 0;
                    clipPolyShape.color = currentColor;
                    clipPolyShape.vertices = clippedPoly;
                    shapes.push_back(clipPolyShape);

                    cout << "Clipped polygon in rectangular window.\n";
                } else {
                    cout << "Clipped polygon is empty (fully outside).\n";
                }
            } else {
                cout << "Not enough vertices for polygon clipping.\n";
            }
            currentPolygon.clear();
        }
        else if (currentShape == ID_CLIP_SQUARE_WINDOW) {
            sqXmin = xc;  sqYmin = yc;
            int dx = xe - xc, dy = ye - yc;
            sqSide = max(abs(dx), abs(dy));
            if (dx < 0) sqXmin = xc - sqSide;
            if (dy < 0) sqYmin = yc - sqSide;
            HPEN hPenS = CreatePen(PS_SOLID, 1, currentColor);
            HGDIOBJ oldPenS = SelectObject(hdc, hPenS);
            Rectangle(hdc, sqXmin, sqYmin, sqXmin + sqSide, sqYmin + sqSide);
            SelectObject(hdc, oldPenS);
            DeleteObject(hPenS);
            shapes.push_back({ "ClipSquareWindow", sqXmin, sqYmin, sqSide, 0, 0, 0, 0, currentColor });
            cout << "Defined square clip window: top-left (" << sqXmin << "," << sqYmin << "), side=" << sqSide << "\n";
        }
        else if (currentShape == ID_CLIP_SQUARE_POINT) {
            bool inside = (xc >= sqXmin && xc <= sqXmin + sqSide && yc >= sqYmin && yc <= sqYmin + sqSide);
            SetPixel(hdc, xc, yc, inside ? RGB(0, 0, 255) : RGB(255, 0, 0));
            shapes.push_back({ "ClipSquarePoint", xc, yc, 0, 0, 0, 0, 0, inside ? RGB(0, 0, 255) : RGB(255, 0, 0) });
            cout << "Square clip point (" << xc << "," << yc << ") is " << (inside ? "inside\n" : "outside\n");
        }
        else if (currentShape == ID_CLIP_SQUARE_LINE) {
            int cx1 = 0, cy1 = 0, cx2 = 0, cy2 = 0;
            if (cohenSutherlandClip(xc, yc, xe, ye, sqXmin, sqYmin, sqXmin + sqSide, sqYmin + sqSide, cx1, cy1, cx2, cy2)) {
                lineBresenham(hdc, cx1, cy1, cx2, cy2, currentColor);
                Shape clipSqLineShape;
                clipSqLineShape.type = "ClipSquareLine";
                clipSqLineShape.x1 = cx1; clipSqLineShape.y1 = cy1;
                clipSqLineShape.x2 = cx2; clipSqLineShape.y2 = cy2;
                clipSqLineShape.a = sqXmin; clipSqLineShape.b = sqYmin;
                clipSqLineShape.c = sqSide; clipSqLineShape.color = currentColor;
                shapes.push_back(clipSqLineShape);
                cout << "Clipped line in square window.\n";
            } else {
                cout << "Line entirely outside square window.\n";
            }
        }
        break;

    case WM_RBUTTONUP:
        if (currentShape == 20) {
            xe2 = LOWORD(lParam);
            ye2 = HIWORD(lParam);
            R  = (int)round(sqrt(pow(xe - xc, 2) + pow(ye - yc, 2)));
            R2 = (int)round(sqrt(pow(xe2 - xc, 2) + pow(ye2 - yc, 2)));
            ellipsePolar(hdc, xc, yc, R, R2, currentColor);
            shapes.push_back({ "EllipsePolar", xc, yc, 0, 0, R, R2, 0, currentColor });
            cout << "Drew Polar ellipse: center (" << xc << "," << yc << "), r1=" << R << ", r2=" << R2 << "\n";
        }
        break;

    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdcPaint = BeginPaint(hwnd, &ps);

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
            else if (s.type == "HermiteSquare") {
                POINT topLeft = { s.x1, s.y1 };
                FillHermiteSquare(hdcPaint, topLeft, s.a, s.color);
            }
            else if (s.type == "FillLines") {
                fillWithLines(hdcPaint, s.x1, s.y1, 0, 0, s.a, s.color, s.b);
            }
            else if (s.type == "FillCircles") {
                fillWithCircles(hdcPaint, s.x1, s.y1, 0, 0, s.a, s.color, s.b);
            }
            else if (s.type == "ClipWindow") {
                circlePolar(hdcPaint, s.x1, s.y1, s.a, s.color);
            }
            else if (s.type == "ClipLine") {
                lineBresenham(hdcPaint, s.x1, s.y1, s.x2, s.y2, s.color);
            }
            else if (s.type == "ClipPoint") {
                SetPixel(hdcPaint, s.x1, s.y1, s.color);
            }
            else if (s.type == "FloodRec") {
                COLORREF targetColor = GetPixel(hdcPaint, s.x1, s.y1);
                floodFillRec(hdcPaint, s.x1, s.y1, targetColor, s.color);
            }
            else if (s.type == "FloodStack") {
                COLORREF targetColor = GetPixel(hdcPaint, s.x1, s.y1);
                floodFillStack(hdcPaint, s.x1, s.y1, targetColor, s.color);
            }
            else if (s.type == "FillBezierRect") {
                POINT topLeft     = { s.x1, s.y1 };
                POINT bottomRight = { s.x2, s.y2 };
                FillBezierRectangle(hdcPaint, topLeft, bottomRight, s.color);
            }
            else if (s.type == "ConvexPoly") {
                fillConvexPolygon(hdcPaint, s.vertices, s.color);
            }
            else if (s.type == "NonConvexPoly") {
                fillNonConvexPolygon(hdcPaint, s.vertices, s.color);
            }
            else if (s.type == "CardinalSpline") {
                drawCardinalSpline(hdcPaint, s.vertices, s.color);
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
                if (!s.vertices.empty()) {
                    HPEN hPenC = CreatePen(PS_SOLID, 1, s.color);
                    HGDIOBJ oldPenC = SelectObject(hdcPaint, hPenC);
                    MoveToEx(hdcPaint, s.vertices[0].x, s.vertices[0].y, NULL);
                    for (size_t i = 1; i < s.vertices.size(); ++i) {
                        LineTo(hdcPaint, s.vertices[i].x, s.vertices[i].y);
                    }
                    LineTo(hdcPaint, s.vertices[0].x, s.vertices[0].y);
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
        }

        EndPaint(hwnd, &ps);
        break;
    }

    case WM_CREATE:
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
// WinMain: Entry Point
//-------------------------------------------------------------------------------------------------

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    ShowWindow(GetConsoleWindow(), SW_SHOW);

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

    if (!RegisterClassEx(&wincl))
        return 0;

    HWND hwnd = CreateWindowEx(
        0,
        _T("CodeBlocksWindowsApp"),
        _T("Computer Graphics Project 2025"),
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 544, 375,
        HWND_DESKTOP, NULL, hInstance, NULL
    );

    ShowWindow(hwnd, nCmdShow);

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

    HMENU hFileMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFileMenu, "File");
    AppendMenu(hFileMenu, MF_STRING, 1, "Save");
    AppendMenu(hFileMenu, MF_STRING, 2, "Load");
    AppendMenu(hFileMenu, MF_STRING, 3, "Clear");

    HMENU hLineMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hLineMenu, "Line");
    AppendMenu(hLineMenu, MF_STRING, 11, "DDA");
    AppendMenu(hLineMenu, MF_STRING, 12, "Midpoint");
    AppendMenu(hLineMenu, MF_STRING, 13, "Parametric");

    HMENU hCircleMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCircleMenu, "Circle");
    AppendMenu(hCircleMenu, MF_STRING, 14, "Direct");
    AppendMenu(hCircleMenu, MF_STRING, 15, "Polar");
    AppendMenu(hCircleMenu, MF_STRING, 16, "Iterative Polar");
    AppendMenu(hCircleMenu, MF_STRING, 17, "Midpoint");
    AppendMenu(hCircleMenu, MF_STRING, 18, "Fast Bresenham");

    HMENU hEllipseMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hEllipseMenu, "Ellipse");
    AppendMenu(hEllipseMenu, MF_STRING, 19, "Direct");
    AppendMenu(hEllipseMenu, MF_STRING, 20, "Polar");

    HMENU hCurveMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCurveMenu, "Curve");
    AppendMenu(hCurveMenu, MF_STRING, 60,             "Hermite Square");
    AppendMenu(hCurveMenu, MF_STRING, ID_CARDINAL_SPLINE, "Cardinal Spline");

    HMENU hColorMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hColorMenu, "Color");
    AppendMenu(hColorMenu, MF_STRING, 41, "Red");
    AppendMenu(hColorMenu, MF_STRING, 42, "Black");
    AppendMenu(hColorMenu, MF_STRING, 43, "Blue");
    AppendMenu(hColorMenu, MF_STRING, 44, "White");
    AppendMenu(hColorMenu, MF_STRING, 45, "Pink");

    HMENU hClipCircleMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClipCircleMenu, "Circular Clipping");
    AppendMenu(hClipCircleMenu, MF_STRING, 51, "Window");
    AppendMenu(hClipCircleMenu, MF_STRING, 52, "Line");
    AppendMenu(hClipCircleMenu, MF_STRING, 53, "Point");

    HMENU hClipRectMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClipRectMenu, "Rectangular Clipping");
    AppendMenu(hClipRectMenu, MF_STRING, ID_CLIP_RECT_WINDOW, "Window");
    AppendMenu(hClipRectMenu, MF_STRING, ID_CLIP_RECT_POINT,  "Point");
    AppendMenu(hClipRectMenu, MF_STRING, ID_CLIP_RECT_LINE,   "Line");
    AppendMenu(hClipRectMenu, MF_STRING, ID_CLIP_RECT_POLY,   "Polygon");

    HMENU hClipSquareMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClipSquareMenu, "Square Clipping");
    AppendMenu(hClipSquareMenu, MF_STRING, ID_CLIP_SQUARE_WINDOW, "Window");
    AppendMenu(hClipSquareMenu, MF_STRING, ID_CLIP_SQUARE_POINT,  "Point");
    AppendMenu(hClipSquareMenu, MF_STRING, ID_CLIP_SQUARE_LINE,   "Line");

    HMENU hLineFillingMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hLineFillingMenu, "Line Filling");
    AppendMenu(hLineFillingMenu, MF_STRING, 101, "Quadrant 1");
    AppendMenu(hLineFillingMenu, MF_STRING, 102, "Quadrant 2");
    AppendMenu(hLineFillingMenu, MF_STRING, 103, "Quadrant 3");
    AppendMenu(hLineFillingMenu, MF_STRING, 104, "Quadrant 4");

    HMENU hCircleFillingMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCircleFillingMenu, "Circle Filling");
    AppendMenu(hCircleFillingMenu, MF_STRING, 201, "Quadrant 1");
    AppendMenu(hCircleFillingMenu, MF_STRING, 202, "Quadrant 2");
    AppendMenu(hCircleFillingMenu, MF_STRING, 203, "Quadrant 3");
    AppendMenu(hCircleFillingMenu, MF_STRING, 204, "Quadrant 4");

    HMENU hFloodMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFloodMenu, "Flood Fill");
    AppendMenu(hFloodMenu, MF_STRING, ID_FLOOD_REC,   "Recursive");
    AppendMenu(hFloodMenu, MF_STRING, ID_FLOOD_STACK, "Non-Recursive");

    HMENU hBezierMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hBezierMenu, "Fill Rectangle (Bezier)");
    AppendMenu(hBezierMenu, MF_STRING, ID_FILL_BEZIERRECT, "Horizontal");

    HMENU hPolyMenu = CreateMenu();
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hPolyMenu, "Fill Polygon");
    AppendMenu(hPolyMenu, MF_STRING, ID_FILL_CONVEX_POLY,    "Convex");
    AppendMenu(hPolyMenu, MF_STRING, ID_FILL_NONCONVEX_POLY, "Non-Convex");

    SetMenu(hwnd, hMenu);
}

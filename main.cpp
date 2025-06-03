#if defined(UNICODE) && !defined(_UNICODE)
#define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif

#include <tchar.h>
#include <windows.h>
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

double LERP(double earlier, double later, double t)
{
        return ((1 - t) * earlier) + (t * later);
}
double quadraticbezier(double x, double y, double z, double t)
{
        double a = LERP(x, y, t);
        double b = LERP(y, z, t);
        double c = LERP(a, b, t);
        return c;
}
void lineBresenham(HDC hdc, int x0, int y0, int x1, int y1, COLORREF color)
{
        //In this episode of "Yousef" el-Halawany explaining things, I explain Bresenham's algorithm.
        //Basically, we define line error function.
        //Line error(x,y) = dx(y-y0) - dy(x-x0)
        //We multiply it by two to avoid using floats later.
        //Line error = 2*dx(y-y0) - 2*dy(x-x0)
        //Line error increases by  2*dx when we increase y by 1.
        //Line error increases by -2*dy when we increase x by 1.
        //Line error increases by (2*dx)-(2*dy) when we increase x by 1 and increase y by 1.
        //Assume abs(dx)>abs(dy) and that dx>=0 and that dy>=0;
        //Initial ACTUAL error will be 0, BUT we will instead use the error of the next midpoint.
        //Initial next midpoint is (x0+1,y0+0.5), so initially the error of the next midpoint will be -(2*dy)+(dx)
        //We ALWAYS x++ because abs(dx)>abs(dy).
        //When the error of the next midpoint is less than 0, the actual point is above the next midpoint so you y++ as well.
        //Depending on the error of the next midpoint, you either {x++;} or {x++;y++;}
        //Whenever you x++, midpointError +=-2*dy
        //You get the point?
        //Also pay attention to the direction.

        int dx, dy, x, y;

        int d;  //midPointError
        int d1; //delta of midPointError in the always case
        int d2; //delta of midPointError in the less than case
        int direction=1;

        if (abs(x1-x0) >= abs(y1 - y0))
        {
                //Case: X always increases.
                //initial value of midPointError is -(2*dy)+(dx)
                if (x0 > x1)
                {
                        swap(x0, x1);
                        swap(y0, y1);
                }
                x = x0;
                y = y0;
                dx = x1 - x0;
                dy = y1 - y0;
                if (dy < 0)
                {
                        direction = -1;
                        dy *= direction;
                }
                d = -(2 * dy) + (dx);
                d1 = -2 * dy;
                d2 =  2 * dx;
                while(x<=x1)
                {
                        SetPixel(hdc, x, y, color);
                        if (d<0)
                        {
                                d += d2;
                                y+=direction;
                        }
                        d += d1;
                        x++;
                }
        }
        else
        {
                //Case: Y always increases.
                //initial value of midPointError is -(dy)+(2*dx)
                if (y0 > y1)
                {
                        swap(x0, x1);
                        swap(y0, y1);
                }
                x = x0;
                y = y0;
                dx = x1 - x0;
                dy = y1 - y0;
                if (dx<0)
                {
                        direction = -1;
                        dx *= direction;
                }
                d = -(dy)+(2 * dx);
                d1 =  2 * dx;
                d2 = -2 * dy;

                while (y <= y1)
                {
                        SetPixel(hdc, x, y, color);
                        if (d > 0)
                        {
                                d += d2;
                                x+=direction;
                        }
                        d += d1;
                        y++;
                }
        }


}
double cubicbezier(double w, double x, double y, double z, double t)
{
        double a = LERP(w, x, t);
        double b = LERP(x, y, t);
        double c = LERP(y, z, t);
        double d = LERP(a, b, t);
        double e = LERP(b, c, t);
        double f = LERP(d, e, t);
        return f;
}

/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure(HWND, UINT, WPARAM, LPARAM);


void menu(HWND);
HMENU hMenu;

/*  Make the class name into a global variable  */
TCHAR szClassName[] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)

{
        HWND hwnd;               /* This is the handle for our window */
        MSG messages;            /* Here messages to the application are saved */
        WNDCLASSEX wincl;        /* Data structure for the windowclass */

        /* The Window structure */
        wincl.hInstance = hInstance;
        wincl.lpszClassName = szClassName;
        wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
        wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
        wincl.cbSize = sizeof(WNDCLASSEX);

        /* Use default icon and mouse-pointer */
        wincl.hIcon = LoadIcon(NULL, IDI_APPLICATION);
        wincl.hIconSm = LoadIcon(NULL, IDI_APPLICATION);
        wincl.hCursor = LoadCursor(NULL, IDC_CROSS);
        wincl.lpszMenuName = NULL;                 /* No menu */
        wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
        wincl.cbWndExtra = 0;                      /* structure or the window instance */
        /* Use Windows's default colour as the background of the window */
        wincl.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);


        /* Register the window class, and if it fails quit the program */
        if (!RegisterClassEx(&wincl))
                return 0;

        /* The class is registered, let's create the program*/
        hwnd = CreateWindowEx(
                0,                   /* Extended possibilites for variation */
                szClassName,         /* Classname */
                _T("Computer Graphics Project 2025"),       /* Title Text */
                WS_OVERLAPPEDWINDOW, /* default window */
                CW_USEDEFAULT,       /* Windows decides the position */
                CW_USEDEFAULT,       /* where the window ends up on the screen */
                544,                 /* The programs width */
                375,                 /* and height in pixels */
                HWND_DESKTOP,        /* The window is a child-window to desktop */
                NULL,                /* No menu */
                hInstance,       /* Program Instance handler */
                NULL                 /* No Window Creation data */
        );

        /* Make the window visible on the screen */
        ShowWindow(hwnd, nCmdShow);

        /* Run the message loop. It will run until GetMessage() returns 0 */
        while (GetMessage(&messages, NULL, 0, 0))
        {
                /* Translate virtual-key messages into character messages */
                TranslateMessage(&messages);
                /* Send message to WindowProcedure */
                DispatchMessage(&messages);
        }

        /* The program return-value is 0 - The value that PostQuitMessage() gave */
        return messages.wParam;
}


/*  This function is called by the Windows function DispatchMessage()  */


void pointClip(HDC hdc, int x, int y, int R, int xl, int yl, COLORREF color)
{
        double delta_x = pow(xl - x, 2), delta_y = pow(yl - y, 2);
        double d = sqrt(delta_x + delta_y);
        if (d <= R)
                SetPixel(hdc, x, y, RGB(0, 0, 255));
        else
                SetPixel(hdc, x, y, RGB(255, 0, 0));

}

/////////////////////////////////////////////////////////////////
/////////////////////////line Algorithms////////////////////////

//DDALine
void DDALine(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
        // calculate dx & dy
        int dx = x2 - x1;
        int dy = y2 - y1;

        // calculate steps required for generating pixels
        int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

        // calculate increment in x & y for each steps
        float Xinc = dx / (float)steps;
        float Yinc = dy / (float)steps;

        // Put pixel for each step
        float X = x1;
        float Y = y1;
        for (int i = 0; i <= steps; i++) {
                SetPixel(hdc, X, Y, color);  // put pixel at (X,Y)
                X += Xinc;           // increment in x at each step
                Y += Yinc;           // increment in y at each step

        }
}


//parametricline
void parametricline(HDC hdc, int x1, int y1, int x2, int y2, int xc, int yc, int R, COLORREF color)
{
        double step, stepX, stepY, parameter = 0;
        stepX = 1.0 / (double)(abs(x1 - x2));
        stepY = 1.0 / (double)(abs(y1 - y2));
        step = min(stepX, stepY);
        for (parameter; parameter <= 1; parameter += step)
        {
                SetPixel(hdc,round(LERP(x1,x2,parameter)), round(LERP(y1, y2, parameter)), color);
        }
}


/////////////////////////////////////////////////////////////////
/////////////////////////circle Algorithms////////////////////////


void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{

        SetPixel(hdc, xc + a, yc + b, color);
        SetPixel(hdc, xc - a, yc + b, color);
        SetPixel(hdc, xc - a, yc - b, color);
        SetPixel(hdc, xc + a, yc - b, color);
        SetPixel(hdc, xc + b, yc + a, color);
        SetPixel(hdc, xc - b, yc + a, color);
        SetPixel(hdc, xc - b, yc - a, color);
        SetPixel(hdc, xc + b, yc - a, color);

}
void DrawQuadrant(HDC hdc, int xc, int yc, int a, int b, int quadrant, COLORREF color)
{
        switch (quadrant)
        {
        case 0:
                SetPixel(hdc, xc + a, yc - b, color);
                SetPixel(hdc, xc + b, yc - a, color);
                break;
        case 1:
                SetPixel(hdc, xc - a, yc - b, color);
                SetPixel(hdc, xc - b, yc - a, color);
                break;
        case 2:
                SetPixel(hdc, xc - a, yc + b, color);
                SetPixel(hdc, xc - b, yc + a, color);
                break;
        case 3:
                SetPixel(hdc, xc + a, yc + b, color);
                SetPixel(hdc, xc + b, yc + a, color);
                break;
        default:
                break;
        }
}

//Direct
void circleDirectMethod(HDC hdc, int xc, int yc, int R, COLORREF color)
{
        double x = 0, y = R;
        double R2 = R * R;
        while (x < y)
        {
                Draw8Points(hdc, xc, yc, x, y, color);
                x += 0.1;
                y = sqrt((double)(R2 - x * x));

        }
}



//CirclePolar
void CirclePolar(HDC hdc, int xc, int yc, int R, COLORREF color)
{

        int x = R, y = 0;
        double theta = 0, dtheta = 1.0 / R;
        Draw8Points(hdc, xc, yc, x, y, color);
        while (x > y)
        {
                theta += dtheta;
                x = round(R * cos(theta));
                y = round(R * sin(theta));
                Draw8Points(hdc, xc, yc, x, y, color);
        }

}


//CircleIterative
void CircleIterative(HDC hdc, int xc, int yc, int R, COLORREF color)
{
        double x = R, y = 0;
        double x1 = 0;
        double dtheta = 1.0 / R, c = cos(dtheta), s = sin(dtheta);
        Draw8Points(hdc, xc, yc, x, y, color);
        while (x + 20 > y)
        {

                x1 = (x * c) - (y * s);
                y = (x * s) + (y * c);
                x = x1;
                Draw8Points(hdc, xc, yc, round(x), round(y), color);
        }

}



//midpoint
void midpoint(HDC hdc, int xc, int yc, int r, COLORREF color)
{
        int x = 0;
        int y = r;
        double d = 1 - r;
        while (x < y)
        {
                if (d <= 0) {
                        d = d + 2 * x + 3;
                        x++;
                }
                else {
                        d = d + 2 * (x - y) + 5;
                        x++;
                        y--;
                }

                Draw8Points(hdc, xc, yc, x, y, color);
        }
}


//modified midpoint
void CircleFasterBresenham(HDC hdc, int xc, int yc, int R, COLORREF color)
{
        int x = 0, y = R;
        int d = 1 - R;
        int c1 = 3, c2 = 5 - 2 * R;
        Draw8Points(hdc, xc, yc, x, y, color);
        while (x < y)
        {
                if (d < 0)
                {
                        d += c1;
                        c2 += 2;
                }
                else
                {
                        d += c2;
                        c2 += 4;
                        y--;
                }
                c1 += 2;
                x++;
                Draw8Points(hdc, xc, yc, x, y, color);
        }
}
void CircleFasterBresenhamQuadrant(HDC hdc, int xc, int yc, int R, int quadrant, COLORREF color)
{
        int x = 0, y = R;
        int d = 1 - R;
        int c1 = 3, c2 = 5 - 2 * R;
        DrawQuadrant(hdc, xc, yc, x, y, quadrant, color);
        while (x < y)
        {
                if (d < 0)
                {
                        d += c1;
                        c2 += 2;
                }
                else
                {
                        d += c2;
                        c2 += 4;
                        y--;
                }
                c1 += 2;
                x++;
                DrawQuadrant(hdc, xc, yc, x, y, quadrant, color);
        }

}


/////////////////////////////////////////////////////////////////
/////////////////////////Ellipse Algorithms////////////////////////

void Draw4Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{

        SetPixel(hdc, xc + a, yc + b, color);
        SetPixel(hdc, xc - a, yc + b, color);
        SetPixel(hdc, xc - a, yc - b, color);
        SetPixel(hdc, xc + a, yc - b, color);
}

//EllipsePolar
void EllipsePolar(HDC hdc, int xc, int yc, int R1, int R2, COLORREF color)
{

        int x = R1, y = 0;
        double theta = 0, dtheta = 1.0 / R1;
        Draw4Points(hdc, xc, yc, x, y, color);
        while (x + R1 > y)
        {
                theta += dtheta;
                x = round(R1 * cos(theta));
                y = round(R2 * sin(theta));
                Draw4Points(hdc, xc, yc, x, y, color);
        }

}
//Anyway midpoint algorithms are a type of Bresenham algorithms.
//Bresenham curve drawing algorithms are algorithms that rely on error functions.
//Midpoint curve drawing algorithms assume one of the variables always increases while the other increases only some of the time.
//You measure two things: you measure the initial error of the next midpoint and measure how much that error changes when you change the variables.
//In midpoint line drawing algorithms, for example, error = dx(y-y0) - dy(x-x0)
//I already explained that earlier, though.
//What would the Bresenham algorithm of an ellipse be, then?
//Error   = (x^2)(b^2)+(y^2)(a^2)–(a^2)(b^2)
//There is an easy way and a hard way.
//Easy way is to divide it to two parts depending on slope and implement a midpoint algorithm.
//Hard way is probably not what you want us to do...
//You specifically asked for midpoint.
void EllipseBresenham(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
        //Error of next midpoint.
        int errorMid;
        //Change to error of next midpoint on x++
        int chx;
        //Change to error of next midpoint on y++
        int chy;
        //Change to error of next midpoint on x++ y++
        int chxy;

        //dy/dx = (-(B^2)x)/((A^2)y)
        //We only care about the abs of the slope and only care that it equals 1.
        //Therefore we substitute 1 = ((B^2)x)/((A^2)y)
        //Therefore the ending point for both loops is when (B^2)x=(A^2)y.
        int x = -a, y = 0;
        //Region 0: abs(slope) is huge, b*b*x >= a*a*y, x increases only sometimes, y increases always.
        //next midpoint is x+0.5,y+1
        //initial error is (0.5^2)(b^2)+(1^2)(a^2)–(a^2)(b^2)
        //We multiply all error calculations by 4 to avoid floating point operations.
        //initial error is now (b^2)+4(a^2)-4(a^2)(b^2)
        errorMid = (b * b) + (4 * a * a) - (4 * a * a * b * b);
        //change in midpoint error on x++ = error(x+1.5,y+1)-error(x+0.5,y+1)
        //Of course multiply with 4.
        //4((x+1.5)^2)(b^2)+4((y+1)^2)(a^2)–4(a^2)(b^2)-(4((x+0.5)^2)(b^2)+4((y+1)^2)(a^2)–4(a^2)(b^2))
        //=4((x+1.5)^2)(b^2)-(4((x+0.5)^2)(b^2))
        //=(8x+8)(b*b)
        while (abs(b*b*x)>=abs(a*a*y))
        {

        }

        x = 0; y = -b;
        //Region 1: abs(slope) is tiny, b*b*x <= a*a*y, y increases only sometimes, x increases always.
        while (abs(b * b * x) <= abs(a * a * y))
        {
        }

}

//EllipseDirect
void EllipseDirect(HDC hdc, int xc, int yc, int a, int b, COLORREF color) {

        int x = 0, y = b;
        Draw4Points(hdc, xc, yc, x, y, color);
        while (x * (b * b) < ((a * a) * y))
        {
                x++;
                y = round(sqrt(((b * b) * ((a * a) - (x * x))) / (a * a)));
                Draw4Points(hdc, xc, yc, x, y, color);
        }
        y = 0, x = a;
        Draw4Points(hdc, xc, yc, x, y, color);
        while (x * (b * b) > ((a * a) * y))
        {
                y++;
                x = round(sqrt(((a * a) * ((b * b) - (y * y))) / (b * b)));
                Draw4Points(hdc, xc, yc, x, y, color);
        }
}


//////////////////////////////////////////////////////////
//load from file
void load(HWND hWnd, HDC& hdc)
{
        string fileName = "picture.bmp";
        if (fileName == "")
                return;
        HBITMAP hBitmap;
        hBitmap = (HBITMAP)::LoadImage(NULL, fileName.c_str(), IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
        HDC hLocalDC;
        hLocalDC = CreateCompatibleDC(hdc);
        BITMAP qBitmap;
        int iReturn = GetObject(reinterpret_cast<HGDIOBJ>(hBitmap), sizeof(BITMAP), reinterpret_cast<LPVOID>(&qBitmap));
        HBITMAP hOldBmp = (HBITMAP)SelectObject(hLocalDC, hBitmap);
        BOOL qRetBlit = BitBlt(hdc, 0, 0, qBitmap.bmWidth, qBitmap.bmHeight, hLocalDC, 0, 0, SRCCOPY);
        SelectObject(hLocalDC, hOldBmp);
        DeleteDC(hLocalDC);
        DeleteObject(hBitmap);
}
/////////////////////////////
//save to file

bool HDCToFile(const char* FilePath, HDC Context, RECT Area, uint16_t BitsPerPixel)
{
        uint32_t Width = Area.right - Area.left;
        uint32_t Height = Area.bottom - Area.top;

        BITMAPINFO Info;
        BITMAPFILEHEADER Header;
        memset(&Info, 0, sizeof(Info));
        memset(&Header, 0, sizeof(Header));
        Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        Info.bmiHeader.biWidth = Width;
        Info.bmiHeader.biHeight = Height;
        Info.bmiHeader.biPlanes = 1;
        Info.bmiHeader.biBitCount = BitsPerPixel;
        Info.bmiHeader.biCompression = BI_RGB;
        Info.bmiHeader.biSizeImage = Width * Height * (BitsPerPixel > 24 ? 4 : 3);
        Header.bfType = 0x4D42;
        Header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);


        char* Pixels = NULL;
        HDC MemDC = CreateCompatibleDC(Context);
        HBITMAP Section = CreateDIBSection(Context, &Info, DIB_RGB_COLORS, (void**)&Pixels, 0, 0);
        DeleteObject(SelectObject(MemDC, Section));
        BitBlt(MemDC, 0, 0, Width, Height, Context, Area.left, Area.top, SRCCOPY);
        DeleteDC(MemDC);

        std::fstream hFile(FilePath, std::ios::out | std::ios::binary);
        if (hFile.is_open())
        {
                hFile.write((char*)&Header, sizeof(Header));
                hFile.write((char*)&Info.bmiHeader, sizeof(Info.bmiHeader));
                hFile.write(Pixels, (((BitsPerPixel * Width + 31) & ~31) / 8) * Height);
                hFile.close();
                DeleteObject(Section);
                return true;
        }

        DeleteObject(Section);
        return false;
}

////////////////////////////////////////////////////////////
////Filling////
void FillWithLines(HDC hdc, int xc, int yc, int a, int b, int R, COLORREF color, int q) {



        int x = R, y = 0;
        double theta = 0, dtheta = 1.0/R;
        if (q == 4) {
                lineBresenham(hdc, xc, yc, xc - x, yc - y, color);
        }
        else if (q == 2) {

                lineBresenham(hdc, xc, yc, xc + x, yc + y, color);

        }
        else if (q == 1) {

                lineBresenham(hdc, xc, yc, xc + x, yc - y, color);

        }
        else if (q == 3) {
                lineBresenham(hdc, xc, yc, xc - x, yc + y, color);
        }
        while (x * R > y)
        {
                theta += dtheta;
                x = round(R * cos(theta));
                y = round(R * sin(theta));
                if (q == 4) {
                        lineBresenham(hdc, xc, yc, xc - x, yc - y, color);
                }
                else if (q == 2) {

                        lineBresenham(hdc, xc, yc, xc + x, yc + y, color);

                }
                else if (q == 1) {

                        lineBresenham(hdc, xc, yc, xc + x, yc - y, color);

                }
                else if (q == 3) {

                        lineBresenham(hdc, xc, yc, xc - x, yc + y, color);

                }
        }
}
void FillWithCircles(HDC hdc, int xc, int yc, int a, int b, int R, COLORREF color, int quadrant)
{
        for (int i = 0; i < R; i++)
        {
                CircleFasterBresenhamQuadrant(hdc, xc, yc, i, quadrant, color);
        }
}

////////////////////////////
int xc, yc;
int xe, ye;
int xe2, ye2;
int R, R2;
int windowX, windowY, windowR;



LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
        static int shape;
        static int color = RGB(0, 0, 0);


        HDC hdc = GetDC(hwnd);
        switch (message)                  /* handle the messages */
        {

        case WM_COMMAND:

                cout << wParam << endl;
                if (wParam == 41) {
                        color = RGB(255, 0, 0);

                }
                else if (wParam == 42) {

                        color = RGB(0, 0, 0);


                }
                else if (wParam == 43) {

                        color = RGB(0, 0, 255);


                }
                else if (wParam == 44) {

                        color = RGB(255, 255, 255);

                }
                else if (wParam == 45) {

                        color = RGB(255, 63, 127);

                }

                else if (wParam == 1) {

                        RECT rect;
                        if (GetWindowRect(hwnd, &rect)) {

                                rect.top += 8;
                                rect.left += 8;
                                HDCToFile("picture.bmp", hdc, rect, 24);
                                ReleaseDC(hwnd, hdc);


                        }
                }

                else if (wParam == 2) {

                        load(hwnd, hdc);

                }
                else if (wParam == 3) {

                        InvalidateRect(hwnd, NULL, TRUE);
                }
                else if (wParam == 101) {

                        FillWithLines(hdc, xc, yc, xe, ye, R, color, 1);

                }
                else if (wParam == 102) {

                        FillWithLines(hdc, xc, yc, xe, ye, R, color, 4);

                }
                else if (wParam == 103) {

                        FillWithLines(hdc, xc, yc, xe, ye, R, color, 3);

                }
                else if (wParam == 104) {

                        FillWithLines(hdc, xc, yc, xe, ye, R, color, 2);

                }
                else if (wParam == 201) {

                        FillWithCircles(hdc, xc, yc, xe, ye, R, color, 0);

                }
                else if (wParam == 202) {

                        FillWithCircles(hdc, xc, yc, xe, ye, R, color, 1);

                }
                else if (wParam == 203) {

                        FillWithCircles(hdc, xc, yc, xe, ye, R, color, 2);

                }
                else if (wParam == 204) {

                        FillWithCircles(hdc, xc, yc, xe, ye, R, color, 3);

                }
                else {

                        shape = wParam;
                }

                break;

        case WM_LBUTTONDOWN:
                xc = LOWORD(lParam);
                yc = HIWORD(lParam);

                if (shape == 51) {
                        windowX = xc;
                        windowY = yc;
                }
                else if (shape == 53) {

                        pointClip(hdc, xc, yc, windowR, windowX, windowY, 0);
                }


                break;

        case WM_RBUTTONDOWN:

                xe = LOWORD(lParam);
                ye = HIWORD(lParam);

                ////////////////////////////////////////
                //////////////line//////////////////


                if (shape == 11) {

                        DDALine(hdc, xc, yc, xe, ye, color);


                }
                else if (shape == 12) {


                        lineBresenham(hdc, xc, yc, xe, ye, color);

                }
                else if (shape == 13) {
                        parametricline(hdc, xc, yc, xe, ye, 0, 0, 0, color);

                }

                ////////////////////////////////////////
                //////////////circle//////////////////

                else if (shape == 14) {

                        R = sqrt(pow((xe - xc), 2) + pow((ye - yc), 2));
                        circleDirectMethod(hdc, xc, yc, R, color);

                }
                else if (shape == 15) {

                        R = sqrt(pow((xe - xc), 2) + pow((ye - yc), 2));
                        CirclePolar(hdc, xc, yc, R, color);
                }

                else if (shape == 16) {

                        R = sqrt(pow((xe - xc), 2) + pow((ye - yc), 2));
                        CircleIterative(hdc, xc, yc, R, color);
                }
                else if (shape == 17) {

                        R = sqrt(pow((xe - xc), 2) + pow((ye - yc), 2));
                        midpoint(hdc, xc, yc, R, color);

                }
                else if (shape == 18) {

                        R = sqrt(pow((xe - xc), 2) + pow((ye - yc), 2));
                        CircleFasterBresenham(hdc, xc, yc, R, color);

                }

                ////////////////////////////////////////
                //////////////Direct Ellipse//////////////////
                else if (shape == 19) {

                        EllipseDirect(hdc, xc, yc, xe, ye, color);

                }

                ////////////////////////////////////////
                //////////////Clipping//////////////////
                else if (shape == 51) {
                        R = sqrt(pow((xe - xc), 2) + pow((ye - yc), 2));
                        windowR = R;
                        CirclePolar(hdc, xc, yc, R, color);
                }
                else if (shape == 52) {

                        parametricline(hdc, xc, yc, xe, ye, windowX, windowY, windowR, color);

                }

                break;

        case WM_RBUTTONUP:  //second point

                ////////////////////////////////////////
                //////////////Polar Ellipse//////////////////

                if (shape == 20) {

                        xe2 = LOWORD(lParam);
                        ye2 = HIWORD(lParam);
                        //first radius
                        R = sqrt(pow((xe - xc), 2) + pow((ye - yc), 2));

                        //second radius
                        R2 = (int)sqrt(pow(xe2 - xc, 2.0) + pow(ye2 - yc, 2.0));
                        //draw Ellipse with Polar algo
                        EllipsePolar(hdc, xc, yc, R, R2, color);

                }
                break;



        case WM_CREATE:
                menu(hwnd);
                break;

        case WM_DESTROY:
                PostQuitMessage(0);       /* send a WM_QUIT to the message queue */
                break;

        default:                      /* for messages that we don't deal with */
                return DefWindowProc(hwnd, message, wParam, lParam);
        }

        return 0;
}




void menu(HWND hwnd)
{

        hMenu = CreateMenu();
        HMENU hFileMenu = CreateMenu();
        //file menu start from 1
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hFileMenu, "File");
        AppendMenu(hFileMenu, MF_STRING, 1, "Save");
        AppendMenu(hFileMenu, MF_STRING, 2, "Load");
        AppendMenu(hFileMenu, MF_STRING, 3, "Clear");


        //line menu
        HMENU hLineleMenu = CreateMenu();
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hLineleMenu, "Line");
        AppendMenu(hLineleMenu, MF_STRING, 11, "DDA");
        AppendMenu(hLineleMenu, MF_STRING, 12, "Midpoint");
        AppendMenu(hLineleMenu, MF_STRING, 13, "Parametric");

        //circle menu
        HMENU hCircleMenu = CreateMenu();
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCircleMenu, "Circle");
        AppendMenu(hCircleMenu, MF_STRING, 14, "Direct");
        AppendMenu(hCircleMenu, MF_STRING, 15, "Polar");
        AppendMenu(hCircleMenu, MF_STRING, 16, "Iterative Polar");
        AppendMenu(hCircleMenu, MF_STRING, 17, "Midpoint");
        AppendMenu(hCircleMenu, MF_STRING, 18, "Modified Midpoint");


        //Ellipse menu
        HMENU hEllipse = CreateMenu();
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hEllipse, "Ellipse");
        AppendMenu(hEllipse, MF_STRING, 19, "Direct");
        AppendMenu(hEllipse, MF_STRING, 20, "Polar");



        //color menu
        HMENU hColor = CreateMenu();
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hColor, "Color");
        AppendMenu(hColor, MF_STRING, 41, "Red");
        AppendMenu(hColor, MF_STRING, 42, "Black");
        AppendMenu(hColor, MF_STRING, 43, "Blue");
        AppendMenu(hColor, MF_STRING, 44, "White");
        AppendMenu(hColor, MF_STRING, 45, "Pretty >:3c");


        //point and line clipping menu
        HMENU hClipping = CreateMenu();
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hClipping, "Clipping");
        AppendMenu(hClipping, MF_STRING, 51, "Window");
        AppendMenu(hClipping, MF_STRING, 52, "Line");
        AppendMenu(hClipping, MF_STRING, 53, "Point");

        HMENU hLineFilling = CreateMenu();
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hLineFilling, "Line Filling");
        AppendMenu(hLineFilling, MF_STRING, 101, "1");
        AppendMenu(hLineFilling, MF_STRING, 102, "2");
        AppendMenu(hLineFilling, MF_STRING, 103, "3");
        AppendMenu(hLineFilling, MF_STRING, 104, "4");

        HMENU hCircleFilling = CreateMenu();
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hCircleFilling, "Circle Filling");
        AppendMenu(hCircleFilling, MF_STRING, 201, "1");
        AppendMenu(hCircleFilling, MF_STRING, 202, "2");
        AppendMenu(hCircleFilling, MF_STRING, 203, "3");
        AppendMenu(hCircleFilling, MF_STRING, 204, "4");




        SetMenu(hwnd, hMenu);
}

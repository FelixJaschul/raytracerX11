#ifndef X11_WRAPPER_H
#define X11_WRAPPER_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>

#ifdef __cplusplus
extern "C" {
#endif

// Window Struct
typedef struct xWindow 
{
    Display *display;
    Window window;
    int screen;
    GC gc;
    int width;
    int height;
    int x, y;
    const char *title;
    uint32_t fcolor;
    uint32_t bcolor;

    XImage *image;
    uint32_t *buffer;
} xWindow;

// Initialize the window struct with default values
static inline void xInit(xWindow *w) 
{
    w->display = XOpenDisplay(NULL);
    w->screen = 0;
    w->window = 0;
    w->gc = 0;
    w->width = 800;
    w->height = 600;
    w->x = 100;
    w->y = 100;
    w->title = "change by win.(title) = ...; <- works for everything else aswell";
    w->fcolor = 0x000000;
    w->bcolor = 0xFFFFFF;
    w->image = NULL;
    w->buffer = NULL;
}

// Window Management
static inline bool xCreateWindow(xWindow *w) 
{
    if (!w->display) return false;

    w->window = XCreateSimpleWindow(
        w->display,
        RootWindow(w->display, w->screen),
        w->x, w->y,
        w->width, w->height, 0,
        BlackPixel(w->display, w->screen),
        WhitePixel(w->display, w->screen)
    );

    if (!w->window) return false;

    XStoreName(w->display, w->window, w->title);
    XSelectInput(w->display, w->window, ExposureMask | KeyPressMask | KeyReleaseMask | StructureNotifyMask);

    Atom wmDeleteMessage = XInternAtom(w->display, "WM_DELETE_WINDOW", False);
    XSetWMProtocols(w->display, w->window, &wmDeleteMessage, 1);
    XMapWindow(w->display, w->window);
    XFlush(w->display);

    XEvent ev;
    XWindowEvent(w->display, w->window, ExposureMask, &ev);
    w->gc = DefaultGC(w->display, w->screen);

    // Create Framebuffer
    w->buffer = (uint32_t *)malloc(w->width * w->height * sizeof(uint32_t));
    w->image = XCreateImage(
        w->display,
        DefaultVisual(w->display, w->screen),
        DefaultDepth(w->display, w->screen),
        ZPixmap, 0, (char *)w->buffer,
        w->width, w->height, 32, 0
    );

    return true;
}

static inline void xDestroyWindow(xWindow *w) 
{
    if (!w->display || !w->window) return;
    if (w->image) XDestroyImage(w->image); // This also frees w->buffer
    XDestroyWindow(w->display, w->window);
    XSync(w->display, False);
    w->window = 0;
}

// Drawing
static inline void xDrawPixel(const xWindow *w, const int x, const int y, const uint32_t color)
{
    if (x >= 0 && x < w->width && y >= 0 && y < w->height)
    {
        w->buffer[y * w->width + x] = color;
    }
}

static inline void xUpdateFramebuffer(const xWindow *w)
{
    XPutImage(w->display, w->window, w->gc, w->image, 0, 0, 0, 0, w->width, w->height);
    XFlush(w->display);
}

static inline void xDrawTriangle(
    const xWindow *w,
    const int x1,
    const int y1,
    const int x2,
    const int y2,
    const int x3,
    const int y3,
    const uint32_t color)
{
    XSetForeground(w->display, w->gc, color);
    XPoint pts[3] = {
        { (short)x1, (short)y1 },
        { (short)x2, (short)y2 },
        { (short)x3, (short)y3 }
    };
    XFillPolygon(w->display, w->window, w->gc, pts, 3, Convex, CoordModeOrigin);
    XFlush(w->display);
}

static inline void xDrawRectangle(
    const xWindow *w,
    const int x1,
    const int y1,
    const int width,
    const int height,
    const uint32_t color)
{
    const int x2 = x1 + width;
    const int y2 = y1 + height;
    xDrawTriangle(w, x1, y1, x2, y1, x1, y2, color);
    xDrawTriangle(w, x2, y1, x2, y2, x1, y2, color);
}

// Input
typedef enum xKey 
{
    Unknown = 0,

    Escape, Space, 
    Enter, Tab, 
    Backspace,

    Left, Right, 
    Up, Down,

    A,B,C,D,
    E,F,G,H,
    I,J,K,L,
    M,N,O,P,
    Q,R,S,T,
    U,V,W,X,
    Y,Z,

    Num0,
    Num1,Num2,Num3,
    Num4,Num5,Num6,
    Num7,Num8,Num9,

    LeftShift, RightShift, 
    LeftControl, RightControl,

    LeftAlt, RightAlt, 
    LeftSuper, RightSuper,

    F1, F2,
    KEY_COUNT
} xKey;

#ifndef __cplusplus
static bool s_curr[KEY_COUNT] = {0};
static bool s_prev[KEY_COUNT] = {0};
#else
#include <array>
static std::array<bool, KEY_COUNT> s_curr{};
static std::array<bool, KEY_COUNT> s_prev{};
#endif

static inline void xSetKey(
    const xKey key,
    const bool down)
{
    const size_t idx = (size_t)key;
    if (idx < KEY_COUNT) s_curr[idx] = down;
}

static inline bool xPollEvents(
    Display *display)
{
    const Atom wmDeleteMessage = XInternAtom(display, "WM_DELETE_WINDOW", False);
    bool shouldClose = false;

    while (XPending(display) > 0) 
    {
        XEvent event;
        XNextEvent(display, &event);

        switch (event.type) 
        {
            case ClientMessage:
                if ((Atom)event.xclient.data.l[0] == wmDeleteMessage) shouldClose = true;
                break;

            case KeyPress:
            case KeyRelease:
            {
                bool down = (event.type == KeyPress);
                KeySym sym = XLookupKeysym(&event.xkey, 0);
                xKey k = Unknown;

                switch (sym) 
                {
                    case XK_Escape: k = Escape; break;
                    case XK_space: k = Space; break;
                    case XK_Return: k = Enter; break;
                    case XK_Tab: k = Tab; break;
                    case XK_BackSpace: k = Backspace; break;
                    case XK_Left: k = Left; break;
                    case XK_Right: k = Right; break;
                    case XK_Up: k = Up; break;
                    case XK_Down: k = Down; break;
                    case XK_a: case XK_A: k = A; break;
                    case XK_b: case XK_B: k = B; break;
                    case XK_c: case XK_C: k = C; break;
                    case XK_d: case XK_D: k = D; break;
                    case XK_e: case XK_E: k = E; break;
                    case XK_f: case XK_F: k = F; break;
                    case XK_g: case XK_G: k = G; break;
                    case XK_h: case XK_H: k = H; break;
                    case XK_i: case XK_I: k = I; break;
                    case XK_j: case XK_J: k = J; break;
                    case XK_k: case XK_K: k = K; break;
                    case XK_l: case XK_L: k = L; break;
                    case XK_m: case XK_M: k = M; break;
                    case XK_n: case XK_N: k = N; break;
                    case XK_o: case XK_O: k = O; break;
                    case XK_p: case XK_P: k = P; break;
                    case XK_q: case XK_Q: k = Q; break;
                    case XK_r: case XK_R: k = R; break;
                    case XK_s: case XK_S: k = S; break;
                    case XK_t: case XK_T: k = T; break;
                    case XK_u: case XK_U: k = U; break;
                    case XK_v: case XK_V: k = V; break;
                    case XK_w: case XK_W: k = W; break;
                    case XK_x: case XK_X: k = X; break;
                    case XK_y: case XK_Y: k = Y; break;
                    case XK_z: case XK_Z: k = Z; break;
                    case XK_0: k = Num0; break;
                    case XK_1: k = Num1; break;
                    case XK_2: k = Num2; break;
                    case XK_3: k = Num3; break;
                    case XK_4: k = Num4; break;
                    case XK_5: k = Num5; break;
                    case XK_6: k = Num6; break;
                    case XK_7: k = Num7; break;
                    case XK_8: k = Num8; break;
                    case XK_9: k = Num9; break;
                    case XK_Shift_L: k = LeftShift; break;
                    case XK_Shift_R: k = RightShift; break;
                    case XK_Control_L: k = LeftControl; break;
                    case XK_Control_R: k = RightControl; break;
                    case XK_Alt_L: k = LeftAlt; break;
                    case XK_Alt_R: k = RightAlt; break;
                    case XK_Super_L: k = LeftSuper; break;
                    case XK_Super_R: k = RightSuper; break;
                    case XK_F1: k = F1; break;
                    case XK_F2: k = F2; break;
                    default: k = Unknown; break;
                }

                if (k != Unknown) xSetKey(k, down);
                break;
            }

            default: break;
        }
    }
    return shouldClose;
}

static inline bool xIsKeyDown(
    const xKey key)
{
    const size_t idx = (size_t)key;
    return (idx < KEY_COUNT) ? s_curr[idx] : false;
}

static inline bool xIsKeyPressed(
    const xKey key)
{
    const size_t idx = (size_t)key;
    return (idx < KEY_COUNT) ? (s_curr[idx] && !s_prev[idx]) : false;
}

static inline bool xIsKeyReleased(
    const xKey key)
{
    size_t const idx = (size_t)key;
    return (idx < KEY_COUNT) ? (!s_curr[idx] && s_prev[idx]) : false;
}

static inline void xUpdateInput(void) 
{
#ifndef __cplusplus
    for (size_t i = 0; i < KEY_COUNT; ++i) s_prev[i] = s_curr[i];
#else
    s_prev = s_curr;
#endif
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // X11_WRAPPER_H


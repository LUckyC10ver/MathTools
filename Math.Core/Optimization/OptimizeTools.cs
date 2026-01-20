using System;

namespace MathTools.Core
{
    public partial class Functions
    {
        public static void OPTbracketMinimum(
            ref double a,
            ref double b,
            ref double c,
            out double fa,
            out double fb,
            out double fc,
            Func<double, double> f)
        {
            if (f == null)
            {
                throw new Exception("OPTbracketMinimum: function is required");
            }

            double gold = 1.618034;
            double tiny = 1e-20;
            double glimit = 100.0;
            double cMax = c;

            fa = f(a);
            fb = f(b);

            if (fb > fa)
            {
                Swap(ref a, ref b);
                Swap(ref fb, ref fa);
            }

            c = Math.Min(b + gold * (b - a), cMax);
            fc = f(c);

            while (fb > fc)
            {
                double r = (b - a) * (fb - fc);
                double q = (b - c) * (fb - fa);
                double u = b - ((b - c) * q - (b - a) * r) /
                    (2.0 * Sign(q - r) * Math.Max(Math.Abs(q - r), tiny));
                double ulim = Math.Min(b + glimit * (c - b), cMax);

                if ((b - u) * (u - c) > 0.0)
                {
                    double fu = f(u);
                    if (fu < fc)
                    {
                        a = b;
                        b = u;
                        fa = fb;
                        fb = fu;
                        return;
                    }

                    if (fu > fb)
                    {
                        c = u;
                        fc = fu;
                        return;
                    }

                    u = c + gold * (c - b);
                    fu = f(u);
                }

                if ((c - u) * (u - ulim) > 0.0)
                {
                    double fu = f(u);
                    if (fu < fc)
                    {
                        a = b;
                        b = c;
                        c = u;
                        fa = fb;
                        fb = fc;
                        fc = fu;
                        u = c + gold * (c - b);
                        fu = f(u);
                    }
                }
                else if ((u - ulim) * (ulim - c) >= 0.0)
                {
                    u = ulim;
                    double fu = f(u);
                    a = b;
                    b = c;
                    c = u;
                    fa = fb;
                    fb = fc;
                    fc = fu;
                    continue;
                }
                else
                {
                    u = c + gold * (c - b);
                    double fu = f(u);
                    a = b;
                    b = c;
                    c = u;
                    fa = fb;
                    fb = fc;
                    fc = fu;
                    continue;
                }

                double fu2 = f(u);
                a = b;
                b = c;
                c = u;
                fa = fb;
                fb = fc;
                fc = fu2;
            }
        }

        public static double OPTgoldenSectionSearch(
            double ax,
            double bx,
            double cx,
            Func<double, double> func,
            double tol,
            out double xmin)
        {
            if (func == null)
            {
                throw new Exception("OPTgoldenSectionSearch: function is required");
            }

            double r = 0.61803399;
            double c = 1.0 - r;
            double x0 = ax;
            double x3 = cx;
            double x1;
            double x2;

            if (Math.Abs(cx - bx) > Math.Abs(bx - ax))
            {
                x1 = bx;
                x2 = bx + c * (cx - bx);
            }
            else
            {
                x2 = bx;
                x1 = bx - c * (bx - ax);
            }

            double f1 = func(x1);
            double f2 = func(x2);

            while (Math.Abs(x3 - x0) > tol * (Math.Abs(x1) + Math.Abs(x2)))
            {
                if (f2 < f1)
                {
                    Shift3(ref x0, ref x1, ref x2, r * x1 + c * x3);
                    Shift2(ref f1, ref f2, func(x2));
                }
                else
                {
                    Shift3(ref x3, ref x2, ref x1, r * x2 + c * x0);
                    Shift2(ref f2, ref f1, func(x1));
                }
            }

            if (f1 < f2)
            {
                xmin = x1;
                return f1;
            }

            xmin = x2;
            return f2;
        }

        public static double OPTbrentsMinimumSearch(
            double ax,
            double bx,
            double cx,
            Func<double, double> f,
            double tol,
            out double xmin,
            int ITMAX = 100)
        {
            if (f == null)
            {
                throw new Exception("OPTbrentsMinimumSearch: function is required");
            }

            double cgold = 0.3819660;
            double zeps = 1.0e-10;

            double a = Math.Min(ax, cx);
            double b = Math.Max(ax, cx);
            double x = bx;
            double w = bx;
            double v = bx;
            double fw = f(x);
            double fv = fw;
            double fx = fw;
            double e = 0.0;
            double d = 0.0;

            for (int iter = 1; iter <= ITMAX; iter++)
            {
                double xm = 0.5 * (a + b);
                double tol1 = tol * Math.Abs(x) + zeps;
                double tol2 = 2.0 * tol1;

                if (Math.Abs(x - xm) <= (tol2 - 0.5 * (b - a)))
                {
                    xmin = x;
                    return fx;
                }

                if (Math.Abs(e) > tol1)
                {
                    double r = (x - w) * (fx - fv);
                    double q = (x - v) * (fx - fw);
                    double p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    if (q > 0.0)
                    {
                        p = -p;
                    }
                    q = Math.Abs(q);
                    double etemp = e;
                    e = d;
                    if (Math.Abs(p) >= Math.Abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
                    {
                        d = cgold * (e = x >= xm ? a - x : b - x);
                    }
                    else
                    {
                        d = p / q;
                        double u = x + d;
                        if (u - a < tol2 || b - u < tol2)
                        {
                            d = Sign(xm - x) * Math.Abs(tol1);
                        }
                    }
                }
                else
                {
                    d = cgold * (e = x >= xm ? a - x : b - x);
                }

                double u2 = Math.Abs(d) >= tol1 ? x + d : x + Sign(d) * Math.Abs(tol1);
                double fu = f(u2);

                if (fu <= fx)
                {
                    if (u2 >= x)
                    {
                        a = x;
                    }
                    else
                    {
                        b = x;
                    }
                    v = w;
                    w = x;
                    x = u2;
                    fv = fw;
                    fw = fx;
                    fx = fu;
                }
                else
                {
                    if (u2 < x)
                    {
                        a = u2;
                    }
                    else
                    {
                        b = u2;
                    }
                    if (fu <= fw || w == x)
                    {
                        v = w;
                        w = u2;
                        fv = fw;
                        fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w)
                    {
                        v = u2;
                        fv = fu;
                    }
                }
            }

            throw new Exception("Too many iterations in OPTbrentsMinimumSearch");
        }

        private static void Shift2(ref double a, ref double b, double c)
        {
            a = b;
            b = c;
        }

        private static void Shift3(ref double a, ref double b, ref double c, double d)
        {
            a = b;
            b = c;
            c = d;
        }

        private static void Swap(ref double a, ref double b)
        {
            double t = a;
            a = b;
            b = t;
        }

        private static double Sign(double value)
        {
            return value >= 0.0 ? 1.0 : -1.0;
        }
    }
}

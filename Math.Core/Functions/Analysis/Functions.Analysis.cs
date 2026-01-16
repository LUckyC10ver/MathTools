using System;

namespace Math.Core
{
    public partial class Functions
    {
        /// <summary>
        /// 用 Neville 算法进行多项式插值。
        /// </summary>
        /// <param name="matrix">每行包含 x、y 的二维数组。</param>
        /// <param name="x">目标 x。</param>
        /// <param name="error">插值误差估计。</param>
        /// <returns>插值结果。</returns>
        public static double polint(double[][] matrix, double x, out double error)
        {
            if (matrix == null || matrix.Length == 0 || matrix[0].Length < 2)
            {
                throw new Exception("matrix is empty or irregular");
            }

            int n = matrix.Length;
            var c = new double[n];
            var d = new double[n];

            int ns = 0;
            double dif = Math.Abs(x - matrix[0][0]);
            for (int i = 0; i < n; i++)
            {
                double dift = Math.Abs(x - matrix[i][0]);
                if (dift < dif)
                {
                    ns = i;
                    dif = dift;
                }

                c[i] = matrix[i][1];
                d[i] = matrix[i][1];
            }

            double result = matrix[ns][1];
            error = 0.0;

            for (int m = 1; m < n; m++)
            {
                for (int i = 0; i < n - m; i++)
                {
                    double ho = matrix[i][0] - x;
                    double hp = matrix[i + m][0] - x;
                    double w = c[i + 1] - d[i];
                    double den = ho - hp;
                    if (den == 0.0)
                    {
                        throw new Exception("two X-Values are identical");
                    }

                    den = w / den;
                    d[i] = hp * den;
                    c[i] = ho * den;
                }

                if (ns > 0)
                {
                    ns--;
                    error = d[ns];
                }
                else
                {
                    error = c[0];
                }

                result += error;
            }

            return result;
        }

        /// <summary>
        /// 梯形法积分的第 n 次细分结果。
        /// </summary>
        public static double integralTrapez(Func<double, double> func, double a, double b, int n, double lastInt)
        {
            n = Math.Max(n, 1);
            if (n == 1)
            {
                return 0.5 * (b - a) * (func(a) + func(b));
            }

            int it = 1;
            for (int j = 1; j < n - 1; j++)
            {
                it *= 2;
            }

            double del = (b - a) / it;
            double x = a + 0.5 * del;
            double sum = 0.0;
            for (int j = 1; j <= it; j++)
            {
                sum += func(x);
                x += del;
            }

            return 0.5 * (lastInt + (b - a) * sum / it);
        }

        /// <summary>
        /// 梯形法积分。
        /// </summary>
        public static double integralTrap(Func<double, double> func, double a, double b)
        {
            double result = integralTrapez(func, a, b, 1, 0.0);
            double resultOld = result;
            int i = 1;

            do
            {
                i++;
                resultOld = result;
                result = integralTrapez(func, a, b, i, resultOld);
                if (i > MathCoreInfo.BcFunctionalAlgorithmsMaxIterations)
                {
                    throw new Exception($"calculation aborted after {i} steps  error = {Math.Abs(resultOld - result)}");
                }
            }
            while (Math.Abs(resultOld - result) > Math.Abs(MathCoreInfo.BcFunctionalAlgorithmsEpsilon * result));

            return result;
        }

        /// <summary>
        /// Simpson 积分。
        /// </summary>
        public static double integralSimp(Func<double, double> func, double a, double b)
        {
            double st = integralTrapez(func, a, b, 1, 0.0);
            double s = -1.0e30;
            double olds;
            double oldst;
            int i = 1;

            do
            {
                i++;
                oldst = st;
                olds = s;
                st = integralTrapez(func, a, b, i, st);
                s = (4.0 * st - oldst) / 3.0;
                if (i > MathCoreInfo.BcFunctionalAlgorithmsMaxIterations)
                {
                    throw new Exception($"calculation aborted after {i} steps  error = {Math.Abs(s - olds)}");
                }
            }
            while (Math.Abs(s - olds) > MathCoreInfo.BcFunctionalAlgorithmsEpsilon * Math.Abs(olds));

            return s;
        }

        /// <summary>
        /// Romberg 积分。
        /// </summary>
        public static double integralRomb(Func<double, double> func, double a, double b)
        {
            const int minLoops = 5;
            int maxIterations = MathCoreInfo.BcFunctionalAlgorithmsMaxIterations;
            var h = new double[maxIterations + 1];
            var s = new double[maxIterations + 1];
            var hAndS = new double[minLoops][];
            for (int k = 0; k < minLoops; k++)
            {
                hAndS[k] = new double[2];
            }

            double ss = 0.0;
            double dss = 0.0;

            h[1] = 1.0;
            s[1] = integralTrapez(func, a, b, 1, 1.0);

            int i;
            for (i = 2; i < minLoops; i++)
            {
                s[i] = integralTrapez(func, a, b, i, s[i - 1]);
                h[i] = h[i - 1] / 4.0;
            }

            i--;

            do
            {
                i++;
                if (i > maxIterations)
                {
                    throw new Exception($"calculation aborted after {i} steps  error = {dss}");
                }

                s[i] = integralTrapez(func, a, b, i, s[i - 1]);
                h[i] = h[i - 1] / 4.0;

                for (int k = 0; k < minLoops; k++)
                {
                    hAndS[k][0] = h[i - 5 + k + 1];
                    hAndS[k][1] = s[i - 5 + k + 1];
                }

                ss = polint(hAndS, 0.0, out dss);
            }
            while (Math.Abs(dss) >= MathCoreInfo.BcFunctionalAlgorithmsEpsilon * Math.Abs(ss)
                   && i < maxIterations
                   && Math.Abs(s[i]) > 1e-10);

            return ss;
        }

        /// <summary>
        /// 中心差分一阶导数。
        /// </summary>
        public static double devLinear(Func<double, double> func, double x, double h)
        {
            if (h <= 0)
            {
                throw new Exception("intervall h equal or less than zero");
            }

            return (func(x + h) - func(x - h)) / (2.0 * h);
        }

        /// <summary>
        /// 递减步长的中心差分导数。
        /// </summary>
        public static double devSimply(Func<double, double> func, double x, double h)
        {
            double result = devLinear(func, x, h);
            double oldResult;
            int i = 1;

            do
            {
                i++;
                oldResult = result;
                h /= 2.0;
                result = devLinear(func, x, h);
                if (i > MathCoreInfo.BcFunctionalAlgorithmsMaxIterations)
                {
                    throw new Exception($"calculation aborted after {i} steps  error = {Math.Abs(result - oldResult)}");
                }
            }
            while (Math.Abs(oldResult - result) > Math.Abs(MathCoreInfo.BcFunctionalAlgorithmsEpsilon * result));

            return result;
        }

        /// <summary>
        /// Ridders 方法求导。
        /// </summary>
        public static double devRiddler(Func<double, double> func, double x, double h)
        {
            int maxIterations = MathCoreInfo.BcFunctionalAlgorithmsMaxIterations;
            var a = new double[maxIterations, maxIterations];

            double error = 1.0e30;
            double result = 0.0;

            a[0, 0] = devLinear(func, x, h);

            int i = 0;
            do
            {
                i++;
                h /= 1.4;
                a[0, i] = devLinear(func, x, h);
                double fac = 1.96;
                for (int j = 1; j <= i; j++)
                {
                    a[j, i] = (a[j - 1, i] * fac - a[j - 1, i - 1]) / (fac - 1.0);
                    fac *= 1.96;
                    double errorT = Math.Max(Math.Abs(a[j - 1, i] - a[j - 1, i - 1]), Math.Abs(a[j, i] - a[j - 1, i - 1]));
                    if (errorT <= error)
                    {
                        error = errorT;
                        result = a[j, i];
                    }
                }

                if (i > maxIterations)
                {
                    throw new Exception($"calculation aborted in dev_ridders() after {i} steps  error = {error}");
                }
            }
            while (Math.Abs(a[i, i] - a[i - 1, i - 1]) < error);

            return result;
        }

        /// <summary>
        /// 计算有限差分步长并返回误差状态。
        /// </summary>
        public static int devOptFDInterval(
            ref double h,
            ref double df,
            ref double d2f,
            Func<double, double> f,
            double x,
            double omega = 1.e-08,
            double eta = 1.e-08,
            int K = 10,
            double epsA = 1.e-08,
            double epsM = 1.e-16,
            double big = 1.e+30)
        {
            const double diffMin = 0.001;
            const double diffMax = 0.1;

            const int checkInitial = 0;
            const int increase = 1;
            const int decrease = 2;
            const int computeEstimate = 3;

            double fOfx = f(x);
            double hbar = 2.0 * (eta + Math.Abs(x)) * Math.Sqrt(epsA / (omega + Math.Abs(fOfx)));

            h = 10.0 * hbar;
            double hs = -1.0;
            double hD2f = -1.0;
            double fOfxPh = f(x + h);
            double fOfxMh = f(x - h);
            df = (fOfxPh - fOfx) / h;
            d2f = (fOfxPh - 2 * fOfx + fOfxMh) / h / h;
            double backDiff = (fOfx - fOfxMh) / h;
            double condErrD2f = Math.Abs(d2f) > epsM ? (4.0 * epsA / h / h / Math.Abs(d2f)) : big;
            double condErrFd = Math.Abs(df) > epsM ? (2.0 * epsA / h / Math.Abs(df)) : big;
            double condErrBd = Math.Abs(backDiff) > epsM ? (2.0 * epsA / h / Math.Abs(backDiff)) : big;

            int calcMode = checkInitial;

            for (int k = 0; k < K; k++)
            {
                switch (calcMode)
                {
                    case checkInitial:
                        if (Math.Max(condErrFd, condErrBd) <= diffMax)
                        {
                            hs = h;
                        }

                        if (condErrD2f < diffMin)
                        {
                            calcMode = decrease;
                        }
                        else if (condErrD2f > diffMax)
                        {
                            calcMode = increase;
                        }
                        else
                        {
                            hD2f = h;
                            calcMode = computeEstimate;
                        }
                        continue;

                    case increase:
                        h *= 10.0;
                        fOfxPh = f(x + h);
                        fOfxMh = f(x - h);
                        df = (fOfxPh - fOfx) / h;
                        d2f = (fOfxPh - 2 * fOfx + fOfxMh) / h;
                        backDiff = (fOfx - fOfxMh) / h;
                        condErrD2f = Math.Abs(d2f) > epsM ? (4.0 * epsA / h / h / Math.Abs(d2f)) : big;
                        condErrFd = Math.Abs(df) > epsM ? (2.0 * epsA / h / Math.Abs(df)) : big;
                        condErrBd = Math.Abs(backDiff) > epsM ? (2.0 * epsA / h / Math.Abs(backDiff)) : big;

                        if (hs < 0.0 && Math.Max(condErrFd, condErrBd) <= diffMax)
                        {
                            hs = h;
                        }

                        if (condErrD2f <= diffMax)
                        {
                            calcMode = computeEstimate;
                        }
                        continue;

                    case decrease:
                    {
                        double hold = h;
                        h /= 10.0;
                        fOfxPh = f(x + h);
                        fOfxMh = f(x - h);
                        df = (fOfxPh - fOfx) / h;
                        d2f = (fOfxPh - 2 * fOfx + fOfxMh) / h;
                        backDiff = (fOfx - fOfxMh) / h;
                        condErrD2f = Math.Abs(d2f) > epsM ? (4.0 * epsA / h / h / Math.Abs(d2f)) : big;
                        condErrFd = Math.Abs(df) > epsM ? (2.0 * epsA / h / Math.Abs(df)) : big;
                        condErrBd = Math.Abs(backDiff) > epsM ? (2.0 * epsA / h / Math.Abs(backDiff)) : big;

                        if (condErrD2f > diffMax)
                        {
                            hD2f = hold;
                            calcMode = computeEstimate;
                        }
                        else
                        {
                            if (Math.Max(condErrFd, condErrBd) <= diffMax)
                            {
                                hs = h;
                            }

                            if (diffMin <= condErrD2f && condErrD2f <= diffMax)
                            {
                                hD2f = h;
                                calcMode = computeEstimate;
                            }
                        }

                        continue;
                    }

                    case computeEstimate:
                    {
                        h = 2.0 * Math.Sqrt(epsA / Math.Abs(condErrD2f));
                        fOfxPh = f(x + h);
                        df = (fOfxPh - fOfx) / h;
                        double ef = 0.5 * h * Math.Abs(d2f) + 2.0 * epsA / h;

                        double f0 = f(x - hD2f);
                        double f1 = f(x + hD2f);
                        double cd = (f1 - f0) / (2 * hD2f);
                        double ebar = Math.Abs(df - cd);

                        if (Math.Max(ef, ebar) <= 0.5 * Math.Abs(df))
                        {
                            return 0;
                        }

                        return -1;
                    }
                }
            }

            if (hs < 0.0)
            {
                h = hbar;
                df = 0.0;
                d2f = 0.0;
                return 1;
            }

            if (condErrD2f > diffMax)
            {
                h = hs;
                fOfxPh = f(x + h);
                df = (fOfxPh - fOfx) / h;
                d2f = 0.0;
                return 2;
            }

            return 3;
        }
    }
}

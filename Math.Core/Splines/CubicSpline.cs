using System;

namespace MathTools.Core
{
    /// <summary>
    /// 三次样条。
    /// </summary>
    public class CubicSpline : SplineBase
    {
        private double[] _d2y = Array.Empty<double>();

        /// <summary>
        /// 默认构造函数。
        /// </summary>
        public CubicSpline()
        {
        }

        /// <summary>
        /// 构造三次样条。
        /// </summary>
        public CubicSpline(double[] xValues, double[] yValues, bool sorted = false, double dyLeft = 1e30, double dyRight = 1e30)
        {
            Init(xValues, yValues, sorted, dyLeft, dyRight);
        }

        /// <summary>
        /// 初始化三次样条。
        /// </summary>
        public void Init(double[] xValues, double[] yValues, bool sorted, double dyLeft = 1e30, double dyRight = 1e30)
        {
            Init(xValues, yValues, sorted);
            int n = _xValues.Length;
            _d2y = new double[n];
            var u = new double[n - 1];

            if (dyLeft > 0.99e30)
            {
                _d2y[0] = 0.0;
                u[0] = 0.0;
            }
            else
            {
                _d2y[0] = -0.5;
                u[0] = 3.0 / (_xValues[1] - _xValues[0]) *
                       ((_yValues[1] - _yValues[0]) / (_xValues[1] - _xValues[0]) - dyLeft);
            }

            for (int i = 1; i < n - 1; i++)
            {
                double sig = (_xValues[i] - _xValues[i - 1]) / (_xValues[i + 1] - _xValues[i - 1]);
                double p = sig * _d2y[i - 1] + 2.0;
                _d2y[i] = (sig - 1.0) / p;
                u[i] = (_yValues[i + 1] - _yValues[i]) / (_xValues[i + 1] - _xValues[i])
                       - (_yValues[i] - _yValues[i - 1]) / (_xValues[i] - _xValues[i - 1]);
                u[i] = (6.0 * u[i] / (_xValues[i + 1] - _xValues[i - 1]) - sig * u[i - 1]) / p;
            }

            double qn;
            double un;
            if (dyRight > 0.99e30)
            {
                qn = 0.0;
                un = 0.0;
            }
            else
            {
                qn = 0.5;
                un = 3.0 / (_xValues[n - 1] - _xValues[n - 2])
                     * (dyRight - (_yValues[n - 1] - _yValues[n - 2]) / (_xValues[n - 1] - _xValues[n - 2]));
            }

            _d2y[n - 1] = (un - qn * u[n - 2]) / (qn * _d2y[n - 2] + 1.0);
            for (int k = n - 2; k >= 0; k--)
            {
                _d2y[k] = _d2y[k] * _d2y[k + 1] + u[k];
            }
        }

        /// <summary>
        /// 计算样条值。
        /// </summary>
        public double Evaluate(double x)
        {
            int n = _xValues.Length;
            int klo = 0;
            int khi = n - 1;
            while (khi - klo > 1)
            {
                int k = (khi + klo) / 2;
                if (_xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            double h = _xValues[khi] - _xValues[klo];
            if (h > _epsilon)
            {
                double a = (_xValues[khi] - x) / h;
                double b = (x - _xValues[klo]) / h;
                double res = a * _yValues[klo] + b * _yValues[khi];
                res += ((a * a * a - a) * _d2y[klo] + (b * b * b - b) * _d2y[khi]) * h * h / 6.0;
                return res;
            }

            return (_yValues[klo] + _yValues[khi]) / 2.0;
        }

        /// <summary>
        /// 计算一阶导数。
        /// </summary>
        public double GetDerivative(double x)
        {
            int n = _xValues.Length;
            int klo = 0;
            int khi = n - 1;
            while (khi - klo > 1)
            {
                int k = (khi + klo) / 2;
                if (_xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            double h = _xValues[khi] - _xValues[klo];
            if (h > _epsilon)
            {
                double a = (_xValues[khi] - x) / h;
                double b = (x - _xValues[klo]) / h;
                double res = (_yValues[khi] - _yValues[klo]) / h;
                res -= (3 * a * a - 1) / 6.0 * h * _d2y[klo];
                res += (3 * b * b - 1) / 6.0 * h * _d2y[khi];
                return res;
            }

            return (_yValues[klo] + _yValues[khi]) / 2.0;
        }

        /// <summary>
        /// 计算二阶导数。
        /// </summary>
        public double GetSecondDerivative(double x)
        {
            int n = _xValues.Length;
            int klo = 0;
            int khi = n - 1;
            while (khi - klo > 1)
            {
                int k = (khi + klo) / 2;
                if (_xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            double h = _xValues[khi] - _xValues[klo];
            if (h > _epsilon)
            {
                double a = (_xValues[khi] - x) / h;
                double b = (x - _xValues[klo]) / h;
                return a * _d2y[klo] + b * _d2y[khi];
            }

            return (_d2y[klo] + _d2y[khi]) / 2.0;
        }

        /// <summary>
        /// 计算三阶导数。
        /// </summary>
        public double GetThirdDerivative(double x)
        {
            int n = _xValues.Length;
            int klo = 0;
            int khi = n - 1;
            while (khi - klo > 1)
            {
                int k = (khi + klo) / 2;
                if (_xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            double h = _xValues[khi] - _xValues[klo];
            if (h > _epsilon)
            {
                return (_d2y[khi] - _d2y[klo]) / h;
            }

            return 0.0;
        }
    }
}

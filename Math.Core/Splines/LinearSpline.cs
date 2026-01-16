using System;

namespace MathTools.Core
{
    /// <summary>
    /// 线性样条。
    /// </summary>
    public class LinearSpline : SplineBase
    {
        /// <summary>
        /// 默认构造函数。
        /// </summary>
        public LinearSpline()
        {
        }

        /// <summary>
        /// 构造线性样条。
        /// </summary>
        public LinearSpline(double[] points, double[] values, bool sorted = false)
        {
            Init(points, values, sorted);
        }

        /// <summary>
        /// 计算样条值。
        /// </summary>
        public double Evaluate(double x)
        {
            if (_xValues.Length == 0 || _yValues.Length == 0)
            {
                throw new Exception("spline has no data");
            }

            int n = IndexFromValue(x);
            if (n + 1 >= _xValues.Length)
            {
                n = _xValues.Length - 2;
            }

            if (n + 1 < _xValues.Length)
            {
                double h = _xValues[n + 1] - _xValues[n];
                if (h > _epsilon)
                {
                    return _yValues[n]
                           + (_yValues[n + 1] - _yValues[n]) * (x - _xValues[n]) / h;
                }

                return (_yValues[n] + _yValues[n + 1]) / 2.0;
            }

            throw new Exception($"argument {x} out of range ({XMin()} - {XMax()})");
        }
    }
}

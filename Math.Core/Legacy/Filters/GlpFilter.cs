using System;

namespace Math.Core.Legacy.Filters
{
    /// <summary>
    /// Gaussian-like低通滤波器：按点的局部权重进行加权平均。
    /// </summary>
    public sealed class GlpFilter
    {
        private readonly double[] _points;
        private readonly double[] _values;
        private readonly double[] _k;

        /// <summary>
        /// 使用控制点初始化。
        /// </summary>
        public GlpFilter(double[] points, double[] values, bool sorted = false)
        {
            if (points == null || values == null || points.Length != values.Length)
            {
                throw new ArgumentException("points and values must have same length");
            }

            _points = (double[])points.Clone();
            _values = (double[])values.Clone();

            if (!sorted)
            {
                Functions.sortPoints(ref _points, ref _values);
            }

            _k = new double[_points.Length];
            InitHalfWidth();
        }

        /// <summary>
        /// 滤波输出。
        /// </summary>
        public double Evaluate(double x)
        {
            double sumWeight = 0.0;
            double sumY = 0.0;
            for (int i = 0; i < _points.Length; i++)
            {
                double help = Math.Exp(-_k[i] * Functions.sqr(_points[i] - x));
                sumWeight += help;
                sumY += help * _values[i];
            }

            return sumY / sumWeight;
        }

        /// <summary>
        /// 计算权重和。
        /// </summary>
        public double WeightOfX(double x)
        {
            double sumWeight = 0.0;
            for (int i = 0; i < _points.Length; i++)
            {
                double help = Math.Exp(-_k[i] * Functions.sqr(_points[i] - x));
                sumWeight += help;
            }

            return sumWeight;
        }

        private void InitHalfWidth()
        {
            if (_points.Length == 0)
            {
                return;
            }

            double averageHalfWidth = 4.0 * (_points[^1] - _points[0]) / _points.Length;
            double kValue = -Math.Log(0.01) / Functions.sqr(averageHalfWidth);
            for (int i = 0; i < _k.Length; i++)
            {
                _k[i] = kValue;
            }
        }
    }
}

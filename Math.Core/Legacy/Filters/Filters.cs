using System;

namespace Math.Core.Legacy.Filters
{
    /// <summary>
    /// Gaussian-like低通滤波器：按点的局部权重进行加权平均。
    /// </summary>
    public sealed class BCGLPFilter
    {
        private readonly double[] _points;
        private readonly double[] _values;
        private readonly double[] _k;

        public BCGLPFilter(double[] points, double[] values, bool sorted = false)
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

        public double WeightOfX(double x)
        {
            double sumWeight = 0.0;
            for (int i = 0; i < _points.Length; i++)
            {
                sumWeight += Math.Exp(-_k[i] * Functions.sqr(_points[i] - x));
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

    /// <summary>
    /// 线性时不变滤波器实现。
    /// </summary>
    public sealed class BCLtiFilter
    {
        private double[] _a = Array.Empty<double>();
        private double[] _b = Array.Empty<double>();
        private double[] _state = Array.Empty<double>();
        private double _out;
        private double _z0;

        public void Init(double[] a, double[] b, double z0)
        {
            if (a == null || b == null || a.Length == 0 || b.Length == 0)
            {
                throw new ArgumentException("coefficient vectors required");
            }

            _a = (double[])a.Clone();
            _b = (double[])b.Clone();
            _z0 = z0;

            int n = Math.Max(_a.Length, _b.Length);
            Array.Resize(ref _a, n);
            Array.Resize(ref _b, n);

            double a0 = _a[0];
            const double eps = 1e-11;
            if (Math.Abs(a0) < eps)
            {
                throw new InvalidOperationException("first coefficient of denominator polynom must not be zero");
            }

            if (Math.Abs(a0 - 1.0) > eps)
            {
                for (int i = 0; i < n; i++)
                {
                    _a[i] /= a0;
                    _b[i] /= a0;
                }
            }

            _state = new double[n - 1];
            for (int i = 0; i < _state.Length; i++)
            {
                _state[i] = _z0;
            }

            _out = _z0;
        }

        public void Reset()
        {
            for (int i = 0; i < _state.Length; i++)
            {
                _state[i] = _z0;
            }

            _out = _z0;
        }

        public double Tick(double input)
        {
            if (_state.Length == 0)
            {
                return input * _b[0];
            }

            _out = input * _b[0] + _state[0];
            int n = _state.Length - 1;
            for (int i = 0; i < n; i++)
            {
                _state[i] = input * _b[i + 1] - _out * _a[i + 1] + _state[i + 1];
            }

            _state[n] = input * _b[n + 1] - _out * _a[n + 1];
            return _out;
        }
    }
}

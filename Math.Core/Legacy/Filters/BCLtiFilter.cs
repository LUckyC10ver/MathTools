using System;

namespace Math.Core.Legacy.Filters
{
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

        /// <summary>
        /// 初始化滤波器。
        /// </summary>
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

        /// <summary>
        /// 重置内部状态。
        /// </summary>
        public void Reset()
        {
            for (int i = 0; i < _state.Length; i++)
            {
                _state[i] = _z0;
            }

            _out = _z0;
        }

        /// <summary>
        /// 计算下一步输出。
        /// </summary>
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
                double next = input * _b[i + 1] - _out * _a[i + 1] + _state[i + 1];
                _state[i] = next;
            }

            _state[n] = input * _b[n + 1] - _out * _a[n + 1];
            return _out;
        }
    }
}

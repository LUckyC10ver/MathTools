using System;

namespace Math.Core.Legacy.Regression
{
    /// <summary>
    /// 一元线性回归：拟合 y = a0 + a1*x。
    /// </summary>
    public sealed class BCTinyLinReg
    {
        private double _sumX;
        private double _sumX2;
        private double _sumY;
        private double _sumYX;
        private int _count;
        private double _a0;
        private double _a1;
        private bool _changed;

        /// <summary>
        /// 添加样本点。
        /// </summary>
        public void AddPoint(double x, double y)
        {
            _sumX += x;
            _sumX2 += x * x;
            _sumY += y;
            _sumYX += y * x;
            _count++;
            _changed = true;
        }

        /// <summary>
        /// 计算回归系数。
        /// </summary>
        public void Calculate()
        {
            if (!_changed || _count == 0)
            {
                return;
            }

            double denom = _count * _sumX2 - _sumX * _sumX;
            if (Math.Abs(denom) > 1e-20)
            {
                _a1 = (_count * _sumYX - _sumY * _sumX) / denom;
                _a0 = (_sumY - _a1 * _sumX) / _count;
                _changed = false;
            }
        }

        /// <summary>
        /// 拟合的截距。
        /// </summary>
        public double A0 => _a0;

        /// <summary>
        /// 拟合的斜率。
        /// </summary>
        public double A1 => _a1;

        /// <summary>
        /// 清空累积统计。
        /// </summary>
        public void Reset()
        {
            _sumX = 0.0;
            _sumX2 = 0.0;
            _sumY = 0.0;
            _sumYX = 0.0;
            _count = 0;
            _a0 = 0.0;
            _a1 = 0.0;
            _changed = false;
        }
    }
}

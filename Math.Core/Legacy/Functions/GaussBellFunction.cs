using System;
using System.Collections.Generic;

namespace Math.Core.Legacy.Functions
{
    /// <summary>
    /// 高斯钟形函数集合，通过权重加权输出。
    /// </summary>
    public abstract class GaussBellFunction
    {
        /// <summary>
        /// 单个 bell 结构。
        /// </summary>
        public sealed class Bell
        {
            private double _x;
            private double _y;
            private double _halfWidth;

            /// <summary>
            /// 创建 bell。
            /// </summary>
            public Bell(double x, double y, double halfWidth)
            {
                _x = x;
                _y = y;
                _halfWidth = Math.Abs(halfWidth);
            }

            /// <summary>
            /// 计算 bell 权重。
            /// </summary>
            public double Weight(double argument)
            {
                const double valueAtHalfWidth = 0.01;
                double distance = Math.Abs(argument - _x);
                if (distance > _halfWidth)
                {
                    return valueAtHalfWidth * _halfWidth / distance;
                }

                return Math.Exp(Math.Log(valueAtHalfWidth) * Functions.sqr(distance / _halfWidth));
            }

            /// <summary>
            /// bell 输出。
            /// </summary>
            public double Evaluate(double argument) => _y * Weight(argument);

            public double X => _x;
            public double Y => _y;
            public double HalfWidth => _halfWidth;

            public void SetX(double value) => _x = value;
            public void SetY(double value) => _y = value;
            public void SetHalfWidth(double value) => _halfWidth = Math.Abs(value);
        }

        private readonly List<Bell> _bells = new();

        /// <summary>
        /// bell 列表。
        /// </summary>
        public IList<Bell> Bells => _bells;

        /// <summary>
        /// 计算函数值。
        /// </summary>
        public double Evaluate(double argument)
        {
            double sum1 = 0.0;
            double sum2 = 0.0;
            foreach (Bell bell in _bells)
            {
                double weight = bell.Weight(argument);
                sum1 += bell.Y * weight;
                sum2 += weight;
            }

            if (Math.Abs(sum2) < 1e-30)
            {
                throw new InvalidOperationException("division by zero");
            }

            return sum1 / sum2;
        }

        /// <summary>
        /// 最小 x。
        /// </summary>
        public double XMin()
        {
            if (_bells.Count == 0)
            {
                throw new InvalidOperationException("size is zero");
            }

            double result = _bells[0].X;
            foreach (Bell bell in _bells)
            {
                result = Math.Min(bell.X, result);
            }

            return result;
        }

        /// <summary>
        /// 最大 x。
        /// </summary>
        public double XMax()
        {
            if (_bells.Count == 0)
            {
                throw new InvalidOperationException("size is zero");
            }

            double result = _bells[0].X;
            foreach (Bell bell in _bells)
            {
                result = Math.Max(bell.X, result);
            }

            return result;
        }

        /// <summary>
        /// 初始化 bell 集合。
        /// </summary>
        public abstract void Init();
    }
}

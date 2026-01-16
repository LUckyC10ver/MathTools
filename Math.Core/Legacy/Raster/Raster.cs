using System;
using System.IO;

namespace Math.Core.Legacy.Raster
{
    /// <summary>
    /// 数值栅格化：将值映射到给定步长的最近网格点。
    /// </summary>
    public sealed class BCRaster
    {
        private double _value;
        private double _raster;

        /// <summary>
        /// 创建 raster。
        /// </summary>
        public BCRaster(double value = 0.0, double raster = 0.0)
        {
            _value = value;
            _raster = raster;
        }

        /// <summary>
        /// 当前值。
        /// </summary>
        public double Value => _value;

        /// <summary>
        /// 设置栅格步长。
        /// </summary>
        public void SetRaster(double raster)
        {
            if (raster < 0.0)
            {
                _raster = 0.0;
                throw new InvalidOperationException("Raster must be positive or zero");
            }

            _raster = raster;
            SetValue(_value);
        }

        /// <summary>
        /// 设置值并进行栅格化。
        /// </summary>
        public void SetValue(double value)
        {
            if (_raster != 0.0)
            {
                double sign = value < 0.0 ? -1.0 : 1.0;
                double scaled = Math.Abs(value / _raster);
                double integer = Math.Floor(scaled);
                double frac = scaled - integer;
                value = frac < 0.5 ? integer * sign * _raster : (integer + 1.0) * sign * _raster;
            }

            _value = value;
        }

        /// <summary>
        /// 输出当前值。
        /// </summary>
        public void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine($"{_value}(+-{0.5 * _raster})");
        }

        /// <summary>
        /// 数值乘法。
        /// </summary>
        public void Multiply(double factor)
        {
            SetValue(_value * factor);
        }

        /// <summary>
        /// 数值相减。
        /// </summary>
        public void Subtract(double value)
        {
            SetValue(_value - value);
        }

        /// <summary>
        /// 数值相加。
        /// </summary>
        public void Add(double value)
        {
            SetValue(_value + value);
        }
    }
}

using System;
using System.IO;

namespace MathTools.Core.Legacy.Raster
{
    /// <summary>
    /// 数值栅格化：将值映射到给定步长的最近网格点。
    /// </summary>
    public sealed class BCRaster
    {
        private double _value;
        private double _raster;

        public BCRaster(double value = 0.0, double raster = 0.0)
        {
            _value = value;
            _raster = raster;
        }

        public double Value => _value;

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

        public void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine($"{_value}(+-{0.5 * _raster})");
        }

        public void Multiply(double factor) => SetValue(_value * factor);
        public void Subtract(double value) => SetValue(_value - value);
        public void Add(double value) => SetValue(_value + value);
    }

    /// <summary>
    /// 布尔值栅格包装。
    /// </summary>
    public sealed class BCboolRaster
    {
        private bool _value;

        public BCboolRaster(bool value) => _value = value;

        public bool Value => _value;

        public void SetValue(bool value) => _value = value;

        public void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine(_value);
        }
    }
}

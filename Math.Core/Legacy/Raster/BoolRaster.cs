using System;
using System.IO;

namespace Math.Core.Legacy.Raster
{
    /// <summary>
    /// 布尔值栅格包装。
    /// </summary>
    public sealed class BoolRaster
    {
        private bool _value;

        /// <summary>
        /// 初始化布尔栅格。
        /// </summary>
        public BoolRaster(bool value)
        {
            _value = value;
        }

        /// <summary>
        /// 当前值。
        /// </summary>
        public bool Value => _value;

        /// <summary>
        /// 设置值。
        /// </summary>
        public void SetValue(bool value)
        {
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

            writer.WriteLine(_value);
        }
    }
}

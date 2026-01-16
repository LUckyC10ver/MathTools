using System.Collections.Generic;

namespace MathTools.Core.Grids
{
    /// <summary>
    /// 多维网格抽象类，提供索引器访问。
    /// </summary>
    public abstract class BCGrid : BCGridBase
    {
        /// <summary>
        /// 通过索引器访问网格值。
        /// </summary>
        public abstract double this[int index] { get; }

        /// <summary>
        /// 直接访问内部值序列。
        /// </summary>
        public abstract IReadOnlyList<double> Values { get; }

        /// <inheritdoc />
        public override double At(int index)
        {
            return this[index];
        }
    }
}

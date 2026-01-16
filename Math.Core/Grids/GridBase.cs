using System;
using System.Collections.Generic;
using System.IO;

namespace MathTools.Core.Grids
{
    /// <summary>
    /// 多维网格基类，定义统一的坐标访问与索引规则。
    /// </summary>
    public abstract class BCGridBase
    {
        /// <summary>
        /// 网格维度。
        /// </summary>
        public abstract int Dimension { get; }

        /// <summary>
        /// 通过一维索引访问网格值。
        /// </summary>
        public abstract double At(int index);

        /// <summary>
        /// 获取指定维度的第 <paramref name="index"/> 个坐标。
        /// </summary>
        public abstract double GetCoordinate(int dimension, int index);

        /// <summary>
        /// 将多维索引映射到一维索引。
        /// </summary>
        public abstract int GetIndex(IReadOnlyList<int> multiIndex);

        /// <summary>
        /// 将一维索引转换为多维索引。
        /// </summary>
        public abstract IReadOnlyList<int> GetMultiindex(int index);

        /// <summary>
        /// 根据坐标获取“小于等于”坐标的最大多维索引。
        /// </summary>
        public abstract IReadOnlyList<int> GetMaximalMultiindex(IReadOnlyList<double> coordinates);

        /// <summary>
        /// 在给定误差范围内匹配坐标，返回对应的多维索引。
        /// </summary>
        public abstract IReadOnlyList<int> GetMultiindex(IReadOnlyList<double> coordinates, double maxError = 2e-2);

        /// <summary>
        /// 获取指定多维索引的网格值。
        /// </summary>
        public abstract double GetPoint(IReadOnlyList<int> multiIndex);

        /// <summary>
        /// 获取指定维度的全部坐标。
        /// </summary>
        public abstract IReadOnlyList<double> GetPoints(int dimension);

        /// <summary>
        /// 网格是否为空。
        /// </summary>
        public abstract bool Empty { get; }

        /// <summary>
        /// 网格最大值。
        /// </summary>
        public abstract double MaxValue();

        /// <summary>
        /// 网格最小值。
        /// </summary>
        public abstract double MinValue();

        /// <summary>
        /// 网格均值（忽略非法值）。
        /// </summary>
        public abstract double MeanValue();

        /// <summary>
        /// 结构化输出。
        /// </summary>
        public abstract void Output(TextWriter writer);

        /// <summary>
        /// 网格点总数。
        /// </summary>
        public abstract int Size { get; }

        /// <summary>
        /// 获取指定维度的点数。
        /// </summary>
        public abstract int SizeOfDimension(int dimension);

        /// <summary>
        /// 指定维度的最大坐标。
        /// </summary>
        public abstract double XMax(int dimension);

        /// <summary>
        /// 指定维度的最小坐标。
        /// </summary>
        public abstract double XMin(int dimension);

        /// <summary>
        /// 从流中读取坐标数据。
        /// </summary>
        public abstract void ReadCoordinates(TextReader reader);

        /// <summary>
        /// 从流中读取网格点与对应值。
        /// </summary>
        public abstract int ReadPoints(TextReader reader);

        /// <summary>
        /// 设置指定维度的坐标。
        /// </summary>
        public abstract void SetPoints(int dimension, IReadOnlyList<double> coord);
    }
}

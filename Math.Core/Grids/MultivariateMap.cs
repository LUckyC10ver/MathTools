using System.IO;

namespace MathTools.Core.Grids
{
    /// <summary>
    /// 多变量映射抽象基类。
    /// </summary>
    public abstract class BCMultivariateMap<TArgs, TResult>
    {
        /// <summary>
        /// 返回类名称。
        /// </summary>
        public abstract string ClassName { get; }

        /// <summary>
        /// 输出对象结构信息。
        /// </summary>
        public abstract void Output(TextWriter writer);

        /// <summary>
        /// 计算映射值。
        /// </summary>
        public abstract TResult Evaluate(TArgs argument);
    }
}

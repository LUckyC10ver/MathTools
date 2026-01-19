using System.IO;

namespace MathTools.Core
{
    /// <summary>
    /// 多变量函数抽象基类。
    /// </summary>
    public abstract class BCMultivariateFunction<TArgs, TResult>
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
        /// 计算函数值。
        /// </summary>
        public abstract TResult Evaluate(TArgs argument);
    }
}

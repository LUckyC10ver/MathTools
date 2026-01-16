using System;
using System.IO;

namespace Math.Core.Polynomials
{
    /// <summary>
    /// 一元函数抽象基类。
    /// </summary>
    public abstract class UnaryFunction
    {
        /// <summary>
        /// 返回类型名称。
        /// </summary>
        public abstract string ClassName { get; }

        /// <summary>
        /// 输出函数的结构信息。
        /// </summary>
        public abstract void Output(TextWriter writer);

        /// <summary>
        /// 获取函数在指定点的值。
        /// </summary>
        public abstract double Invoke(double argument);

        /// <summary>
        /// 获取在指定点的导数值（默认数值微分）。
        /// </summary>
        public virtual double GetDerivative(double argument)
        {
            return Functions.devSimply(Invoke, argument, 0.1);
        }

        /// <summary>
        /// 获取区间积分（默认梯形法）。
        /// </summary>
        public virtual double GetIntegral(double a, double b)
        {
            return Functions.integralTrap(Invoke, a, b);
        }
    }
}

using System;

namespace MathTools.Console
{
    internal static class Program
    {
        private static void Main(string[] args)
        {
            System.Console.WriteLine("MathTools.Console test harness.");
            System.Console.WriteLine($"OPT version: {MathTools.Core.OPT.GetVersion()}");
            System.Console.WriteLine();

            var registry = OptimizeTests.BuildRegistry();
            if (args.Length == 0)
            {
                System.Console.WriteLine("Available tests:");
                foreach (var name in registry.Keys)
                {
                    System.Console.WriteLine($"  {name}");
                }
                System.Console.WriteLine("Use: MathTools.Console.exe <test-name|all>");
                return;
            }

            if (args.Length == 1 && args[0].Equals("all", StringComparison.OrdinalIgnoreCase))
            {
                OptimizeTests.RunAll(registry);
                return;
            }

            if (!registry.TryGetValue(args[0], out var test))
            {
                System.Console.WriteLine($"Unknown test '{args[0]}'.");
                return;
            }

            OptimizeTests.RunSingle(args[0], test);
        }
    }
}

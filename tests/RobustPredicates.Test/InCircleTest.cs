using System;
using System.Globalization;
using System.IO;
using System.Linq;
using Xunit;

namespace RobustPredicates.Test
{
    public class InCircleTest
    {
        private const int NSimpleData = 1000;
        private const int NPscicoData = 1000;

        [Fact]
        public void Fast_ShouldSucceed()
        {
            Assert.True(InCircle.Fast(
                new double[] { 0, -1 },
                new double[] { 0, 1 },
                new double[] { 1, 0 },
                new double[] { 0, 0 }) < 0);

            Assert.True(InCircle.Fast
                (new double[] { 0, -1 },
                new double[] { 0, 1 },
                new double[] { 1, 0 },
                new double[] { -1, 0 }) == 0);

            Assert.True(InCircle.Fast(
                new double[] { 0, -1 },
                new double[] { 0, 1 },
                new double[] { 1, 0, },
                new double[] { 10, 10 }) > 0);
        }

        [Fact]
        public void Fast_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("test_data/simple_data/incricle-poinst2d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
            double[] results =
             File.ReadAllLines("test_data/simple_data/results-incircle.txt")
             .Select(n => n.Split()).SelectMany(x => x)
             .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

            int count = 0;
            for (int i = 0; i < NSimpleData; i += 8)
            {
                Assert.Equal(Math.Sign(results[count++]), Math.Sign(InCircle.Fast(
                    new double[] { numbers[i], numbers[i + 1] },
                    new double[] { numbers[i + 2], numbers[i + 3] },
                    new double[] { numbers[i + 4], numbers[i + 5] },
                    new double[] { numbers[i + 6], numbers[i + 7] })));
            }
        }

        [Fact]
        public void Robust_ShouldSucceed()
        {
            Assert.True(InCircle.Robust(
                new double[] { 0, -1 },
                new double[] { 0, 1 },
                new double[] { 1, 0 },
                new double[] { 0, 0 }) < 0);
            Assert.True(InCircle.Robust(
                new double[] { 0, -1 },
                new double[] { 0, 1 },
                new double[] { 1, 0 },
                new double[] { -1, 0 }) == 0);
            Assert.True(InCircle.Robust(
                new double[] { 0, -1 },
                new double[] { 0, 1 },
                new double[] { 1, 0, },
                new double[] { 10, 10 }) > 0);
        }

        [Fact]
        public void Robust_FromFile_ShouldSucceed()
        {
            double[] points =
              File.ReadAllLines("test_data/simple_data/incricle-poinst2d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
            double[] results =
             File.ReadAllLines("test_data/simple_data/results-incircle.txt")
             .Select(n => n.Split()).SelectMany(x => x)
             .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

            int count = 0;
            for (int i = 0; i < NSimpleData; i += 8)
            {
                Assert.Equal(Math.Sign(results[count++]), Math.Sign(InCircle.Robust(
                    new double[] { points[i], points[i + 1] },
                    new double[] { points[i + 2], points[i + 3] },
                    new double[] { points[i + 4], points[i + 5] },
                    new double[] { points[i + 6], points[i + 7] })));
            }
        }

        [Fact]
        public void Robust_FromFile_Pscico_ShouldSucceed()
        {
            double[] points =
              File.ReadAllLines("test_data/pscico_data/incircle2d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

            for (int i = 0; i < NPscicoData; i += 10)
            {
                Assert.Equal(Math.Sign(points[i + 9]), Math.Sign(InCircle.Robust(
                    new double[] { points[i + 1], points[i + 2] },
                    new double[] { points[i + 3], points[i + 4] },
                    new double[] { points[i + 5], points[i + 6] },
                    new double[] { points[i + 7], points[i + 8] })));
            }
        }
    }
}

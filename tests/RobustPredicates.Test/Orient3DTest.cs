﻿using System;
using System.Globalization;
using System.IO;
using System.Linq;
using Xunit;

namespace RobustPredicates.Test
{
    public class Orient3DTest
    {
        private const int NSimpleData = 1000;
        private const int NPscicoData = 1000;

        [Fact]
        public void Robust_ShouldSucceed()
        {
            Assert.True(Orient3D.Robust
                (new double[] { 0, 0, 0 },
                new double[] { 1, 0, 0 },
                new double[] { 0, 1, 1 },
                new double[] { 0, 0, 1 }) < 0);
            Assert.True(Orient3D.Robust(
                new double[] { 0, 0, 0 },
                new double[] { 0, 1, 0 },
                new double[] { 1, 0, 0 },
                new double[] { 0, 0, 1 }) > 0);
            Assert.True(Orient3D.Robust(
                new double[] { 0, 0, 0 },
                new double[] { 1, 1, 0 },
                new double[] { 0, 1, 1 },
                new double[] { 0, 1, 1 }) == 0);
        }

        [Fact]
        public void Robust_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("test_data/simple_data/orient3d-points3d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
            double[] results =
             File.ReadAllLines("test_data/simple_data/results-orient3d.txt")
             .Select(n => n.Split()).SelectMany(x => x)
             .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

            int count = 0;
            for (int i = 0; i < NSimpleData; i += 12)
            {
                var orient3DResult = Orient3D.Robust(
                  new double[] { numbers[i], numbers[i + 1], numbers[i + 2] },
                  new double[] { numbers[i + 3], numbers[i + 4], numbers[i + 5] },
                  new double[] { numbers[i + 6], numbers[i + 7], numbers[i + 8] },
                  new double[] { numbers[i + 9], numbers[i + 10], numbers[i + 11] });
                Assert.Equal(Math.Sign(results[count++]), Math.Sign(orient3DResult));
            }
        }

        [Fact]
        public void Fast_ShouldSucceed()
        {
            Assert.True(Orient3D.Fast(
                new double[] { 0, 0, 0 },
                new double[] { 1, 0, 0 },
                new double[] { 0, 1, 1 },
                new double[] { 0, 0, 1 }) < 0);
            Assert.True(Orient3D.Fast(
                new double[] { 0, 0, 0 },
                new double[] { 0, 1, 0 },
                new double[] { 1, 0, 0 },
                new double[] { 0, 0, 1 }) > 0);
            Assert.True(Orient3D.Fast(
                new double[] { 0, 0, 0 },
                new double[] { 1, 1, 0 },
                new double[] { 0, 1, 1 },
                new double[] { 0, 1, 1 }) == 0);
        }

        [Fact]
        public void Fast_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("test_data/simple_data/orient3d-points3d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
            double[] results =
             File.ReadAllLines("test_data/simple_data/results-orient3d.txt")
             .Select(n => n.Split()).SelectMany(x => x)
             .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

            int count = 0;
            for (int i = 0; i < NSimpleData; i += 12)
            {
                var orient3DResult = Orient3D.Fast(
                  new double[] { numbers[i], numbers[i + 1], numbers[i + 2] },
                  new double[] { numbers[i + 3], numbers[i + 4], numbers[i + 5] },
                  new double[] { numbers[i + 6], numbers[i + 7], numbers[i + 8] },
                  new double[] { numbers[i + 9], numbers[i + 10], numbers[i + 11] });
                Assert.Equal(Math.Sign(results[count++]), Math.Sign(orient3DResult));
            }
        }

        [Fact]
        public void Robust_FromFile_Pscico_ShouldSucceed()
        {
            double[] points =
              File.ReadAllLines("test_data/pscico_data/orient3d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

            for (int i = 0; i < NPscicoData; i += 14)
            {
                Assert.Equal(Math.Sign(points[i + 13]), Math.Sign(Orient3D.Robust(
                 new double[] { points[i + 1], points[i + 2], points[i + 3] },
                 new double[] { points[i + 4], points[i + 5], points[i + 6] },
                 new double[] { points[i + 7], points[i + 8], points[i + 9] },
                 new double[] { points[i + 10], points[i + 11], points[i + 12] })));
            }
        }
    }
}

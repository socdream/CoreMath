using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace CoreMath
{
    public static class SegmentExtensions
    {
        public static float[] SegmentClosestPoint(this float[] a, float[] b, float[] point)
        {
            return a.SegmentClosestPoint(b, point, out float dist);
        }
        public static float[] SegmentClosestPoint(this float[] a, float[] b, float[] point, out float dist){
            var dir = b.Substract(a);

            dist = (point.Substract(a).VectorDotProduct(dir) / dir.VectorLengthSq()).Clamp01();

            return a.Add(dir.VectorScale(dist));
        }
    }
}

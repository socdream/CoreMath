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
        /// <summary>
        /// Computes the closes point from a point to the segment. These need to be 3 components vectors
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="point"></param>
        /// <param name="dist"></param>
        /// <returns></returns>
        public static float[] SegmentClosestPoint(this float[] a, float[] b, float[] point, out float dist){
            var dir = b.Substract(a);

            dist = (point.Substract(a).VectorDotProduct(dir) / dir.VectorLengthSq()).Clamp01();

            return a.Add(dir.VectorScale(dist));
        }
        public static float SegmentPointDistanceSquared(this float[] a, float[] b, float[] point)
        {
            var n = b.Substract(a);
            var pa = a.Substract(point);

            float c = n.VectorDotProduct(pa);

            // Closest point is a
            if (c > 0.0f)
                return pa.VectorDotProduct(pa);

            var bp = point.Substract(b);

            // Closest point is b
            if (n.VectorDotProduct(bp) > 0.0f)
                return bp.VectorDotProduct(bp);

            // Closest point is between a and b
            var e = pa.Substract(n.VectorScale(c / n.VectorDotProduct(n)));

            return e.VectorDotProduct(e);
        }
    }
}

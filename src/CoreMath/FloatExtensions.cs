using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace CoreMath
{
    public static class FloatExtensions
    {
        public static float ToDegrees(this float value)
        {
            return value * 180f / (float)Math.PI;
        }

        public static float ToRadians(this float value)
        {
            return value * (float)Math.PI / 180f;
        }
        
        public static float Clamp(this float value, float min, float max)
        {
            if (value < min)
                return min;

            if (value > max)
                return max;

            return value;
        }

        public static float Clamp01(this float value)
        {
            if (value < 0)
                return 0;

            if (value > 1f)
                return 1f;

            return value;
        }
    }
}

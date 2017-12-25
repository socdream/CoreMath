using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace CoreMath.Test
{
    [TestClass]
    public class MathTests
    {
        [TestMethod]
        public void TestFloatExtensions()
        {
            var value = 10f;
            var clamped = value.Clamp01();

            Assert.AreEqual(1f, clamped);

            clamped = value.Clamp(0, 5f);

            Assert.AreEqual(5f, clamped);

            Assert.AreEqual(180f, ((float)Math.PI).ToDegrees());

            Assert.AreEqual((float)Math.PI, (180f).ToRadians());
        }

        [TestMethod]
        public void TestAngleBetween()
        {
            var a = new[] { 1f, 0, 0 };
            var b = new[] { 0, 1f, 0 };
            var c = new[] { 1f, 1f, 0 };

            var angle1 = a.AngleBetweenVectors(b);
            var expected1 = (float)Math.PI / 2f;

            var angle2 = a.AngleBetweenVectors(c);
            var expected2 = (float)Math.PI / 4f;

            var epsilon = float.Epsilon;
            Assert.IsTrue(angle1 <= expected1 + epsilon && angle1 >= expected1 - epsilon);
            Assert.IsTrue(angle2 <= expected2 + epsilon && angle2 >= expected2 - epsilon);
        }
    }
}

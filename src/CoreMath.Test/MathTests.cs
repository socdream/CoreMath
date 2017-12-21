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
    }
}

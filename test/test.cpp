#include "gtest/gtest.h"
#include "Complex.hpp"

constexpr Complex<double> z1(1.3, -3.2);
constexpr Complex<double> z2(-7.3, 0.3);

constexpr double value(9.56);
constexpr double precision(1e-7);

#define ASSERT_COMPLEX_DOUBLE_EQ(lhs, rhs){                 \
            ASSERT_NEAR(lhs.real(), rhs.real(), precision); \
            ASSERT_NEAR(lhs.imag(), rhs.imag(), precision); \
}

TEST(Construct, Default) 
{
    Complex<double> z;

    ASSERT_NEAR(0.0, z.real(), precision);
    ASSERT_NEAR(0.0, z.imag(), precision);
}

TEST(Construct, RealPart)
{
    Complex<double> z(1.2);

    ASSERT_NEAR(1.2, z.real(), precision);
    ASSERT_NEAR(0.0, z.imag(), precision);
}

TEST(Construct, ImagPart) 
{
    Complex<double> z(0, -3.1);

    ASSERT_NEAR(0.0, z.real(), precision);
    ASSERT_NEAR(-3.1, z.imag(), precision);
}

TEST(Construct, Full)
{
    ASSERT_NEAR(1.3, z1.real(), precision);
    ASSERT_NEAR(-3.2, z1.imag(), precision);

    ASSERT_NEAR(-7.3, z2.real(), precision);
    ASSERT_NEAR(0.3, z2.imag(), precision);
}

TEST(Construct, Copy) 
{
    Complex<double> z(z1);
    
    ASSERT_NEAR(z.real(), z1.real(), precision);
    ASSERT_NEAR(z.imag(), z1.imag(), precision);
}

TEST(Oprator, AssignCopy)
{
    Complex<double> z;
    z = z1;

    ASSERT_NEAR(z.real(), z1.real(), precision);
    ASSERT_NEAR(z.imag(), z1.imag(), precision);
}

// TEST Binary operator==
TEST(CompareEqual, ComplexType)
{
    EXPECT_EQ(z1 == z1, true);
    EXPECT_EQ(z1 == z2, false);
}

TEST(CompareEqual, PrimetiveType)
{
    Complex<double> z(9.56, 0);
    EXPECT_EQ(z == value, true);
    EXPECT_EQ(value == z, true);

    EXPECT_EQ(z1 == value, false);
    EXPECT_EQ(value == z1, false);
}

// TEST Binary operator!=
TEST(CompareNotEqual, ComplexType)
{
    EXPECT_EQ(z1 != z1, false);
    EXPECT_EQ(z1 != z2, true);
}

TEST(CompareNotEqual, PrimetiveType)
{
    Complex<double> z(9.56, 0);
    EXPECT_EQ(z != value, false);
    EXPECT_EQ(value != z, false);

    EXPECT_EQ(z1 != value, true);
    EXPECT_EQ(value != z1, true);
}

// TEST Unary operator+
TEST(Unary, Plus)
{
    const auto Result1 = +z1;
    const auto Result2 = +z2;

    ASSERT_COMPLEX_DOUBLE_EQ(Result1, z1);
    ASSERT_COMPLEX_DOUBLE_EQ(Result2, z2);
}

// TEST Unary operator-
TEST(Unary, Minus)
{
    const auto Result1 = -z1;
    const auto Result2 = -z2;

    ASSERT_COMPLEX_DOUBLE_EQ(Result1, Complex<double>(-1.3, 3.2));
    ASSERT_COMPLEX_DOUBLE_EQ(Result2, Complex<double>(7.3, -0.3));
}

// TEST Unary operator+=
TEST(AddAssign, ComplexType)
{
    auto Result = z1;
    Result += z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(-6.0, -2.9));
}

TEST(AddAssign, PrimetiveType)
{
    auto Result = z1;
    Result += value;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(10.86, -3.2));
}

// TEST Binary operator+
TEST(Add_binary, ComplexType)
{
    const auto Result = z1 + z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(-6.0, -2.9));
}

TEST(Add_binary, PrimetiveType)
{
    const auto Result1 = z1 + value;    // operator+(complex, primetive_type)
    const auto Result2 = -value + z1;   // operator+(primetive_type, complex)

    ASSERT_COMPLEX_DOUBLE_EQ(Result1, Complex<double>(10.86, -3.2));
    ASSERT_COMPLEX_DOUBLE_EQ(Result2, Complex<double>(-8.26, -3.2));
}

// TEST Unary operator-=
TEST(SubAssign, ComplexType)
{
    auto Result = z1;
    Result -= z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(8.6, -3.5));
}

TEST(SubAssign, PrimetiveType)
{
    auto Result = z1;
    Result -= value;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(-8.26, -3.2));
}

// TEST Binary operator-
TEST(Sub_binary, ComplexType)
{
    const auto Result = z1 - z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(8.6, -3.5));
}

TEST(Sub_binary, PrimetiveType)
{
    const auto Result1 = z1 - value;    // operator-(complex, primetive_type)
    const auto Result2 = value - z1;    // operator-(+primetive_type, complex)
    const auto Result3 = -value - z1;    // operator-(-primetive_type, complex)

    ASSERT_COMPLEX_DOUBLE_EQ(Result1, Complex<double>(-8.26, -3.2));
    ASSERT_COMPLEX_DOUBLE_EQ(Result2, Complex<double>(8.26, 3.2));
    ASSERT_COMPLEX_DOUBLE_EQ(Result3, Complex<double>(-10.86, 3.2));
}

// TEST Unary operator*=
TEST(MpyAssign, ComplexType)
{
    auto Result = z1;
    Result *= z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(-8.53, 23.75));
}

TEST(MpyAssign, PrimetiveType)
{
    auto Result = z1;
    Result *= value;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(12.428, -30.592));
}

// TEST Binary operator*
TEST(Mpy_binary, ComplexType)
{
    const auto Result = z1 * z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(-8.53, 23.75));
}

TEST(Mpy_binary, PrimetiveType)
{
    const auto Result1 = z1 * value;    // operator*(complex, primetive_type)
    const auto Result2 = (-value) * z1; // operator*(-primetive_type, complex)
    ASSERT_COMPLEX_DOUBLE_EQ(Result1, Complex<double>(12.428, -30.592));
    ASSERT_COMPLEX_DOUBLE_EQ(Result2, Complex<double>(-12.428, 30.592));
}

// TEST Unary operator/=
TEST(DivAssign, ComplexType)
{
    auto Result = z1;
    Result /= z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(-0.195766204571, 0.430310977894342));
}

TEST(DivAssign, PrimetiveType)
{
    auto Result = z1;
    Result /= value;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(0.1359832635983, -0.33472803347));
}

// TEST Binary operator/
TEST(Div_binary, ComplexType)
{
    const auto Result = z1 / z2;
    ASSERT_COMPLEX_DOUBLE_EQ(Result, Complex<double>(-0.195766204571, 0.430310977894342));
}

TEST(Div_binary, PrimetiveType)
{
    const auto Result1 = z1 / value;        // operator/(complex, primetive_type)
    const auto Result2 = (-value) / z1;     // operator/(-primetive_type, complex)
    ASSERT_COMPLEX_DOUBLE_EQ(Result1, Complex<double>(0.1359832635983, -0.33472803347));
    ASSERT_COMPLEX_DOUBLE_EQ(Result2, Complex<double>(-0.1359832635983, 0.33472803347));
}

TEST(Function, Magnitude)
{
    ASSERT_NEAR(abs(z1), 3.45398320783411, precision);
    ASSERT_NEAR(abs(z2), 7.30616178304313, precision);
}

TEST(Function, Phase)
{
    ASSERT_NEAR(arg(z1), -1.1849136574, precision);
    ASSERT_NEAR(arg(z2), 3.100519875, precision);
}

TEST(Function, Norm)
{
    ASSERT_NEAR(norm(z1), 11.93, precision);
    ASSERT_NEAR(norm(z2), 53.38, precision);
}

TEST(Function, Conjugation)
{
    ASSERT_COMPLEX_DOUBLE_EQ(conjugation(z1), Complex<double>(1.3, 3.2));
    ASSERT_COMPLEX_DOUBLE_EQ(conjugation(z2), Complex<double>(-7.3, -0.3));
}

TEST(Function, Real)
{
    ASSERT_NEAR(real(z1), 1.3, precision);
    ASSERT_NEAR(real(z2), -7.3, precision);
}

TEST(Function, Imag)
{
    ASSERT_NEAR(imag(z1), -3.2, precision);
    ASSERT_NEAR(imag(z2), 0.3, precision);
}

TEST(Function, Exp)
{
    ASSERT_COMPLEX_DOUBLE_EQ(exp(z1), Complex<double>(-3.66303969413, 0.214192049954));
    ASSERT_COMPLEX_DOUBLE_EQ(exp(z2), Complex<double>(0.000645366841762, 0.000199635358453));
}

TEST(Function, Log)
{
    ASSERT_COMPLEX_DOUBLE_EQ(log(z1), Complex<double>(1.239528179, -1.184913635));
    ASSERT_COMPLEX_DOUBLE_EQ(log(z2), Complex<double>(1.988718033, 3.100519896));
}

TEST(Function, Log10)
{
    ASSERT_COMPLEX_DOUBLE_EQ(log10(z1), Complex<double>(0.5383202434, -0.514601469));
    ASSERT_COMPLEX_DOUBLE_EQ(log10(z2), Complex<double>(0.8636893034, 1.346538663));
}

TEST(Function, Sqrt)
{
    ASSERT_COMPLEX_DOUBLE_EQ(sqrt(z1), Complex<double>(1.54174952697, -1.03778205993));
    ASSERT_COMPLEX_DOUBLE_EQ(sqrt(z2), Complex<double>(0.0555057791727, 2.70242130163));
}

TEST(Function, Nth_root)
{
    ASSERT_THROW(nth_root(z1, 0), std::invalid_argument);
    
    ASSERT_COMPLEX_DOUBLE_EQ(nth_root(z2, 1), z2);
    ASSERT_COMPLEX_DOUBLE_EQ(nth_root(z1, 1), z1);
    ASSERT_COMPLEX_DOUBLE_EQ(nth_root(z2, -2), Complex<double>(0.007597119, -0.369882489));

    ASSERT_COMPLEX_DOUBLE_EQ(nth_root(z1, 9), Complex<double>(0.968393629, 0.615904164));
    ASSERT_COMPLEX_DOUBLE_EQ(nth_root(z2, 8), Complex<double>(1.1871166223, 0.484594426));
}

TEST(Function, Make_real)
{
    ASSERT_COMPLEX_DOUBLE_EQ(make_complex_real(2.38), Complex<double>(2.38, 0.0));
}

TEST(Function, Make_imag)
{
    ASSERT_COMPLEX_DOUBLE_EQ(make_complex_imag(-2.38), Complex<double>(0.0, -2.38));
}

TEST(Function, Make_complex)
{
    ASSERT_COMPLEX_DOUBLE_EQ(make_complex(1.3, -3.2), z1);
}

// TEST pow(complex, int)
TEST(Pow, Complex_to_int)
{
    ASSERT_COMPLEX_DOUBLE_EQ(pow(z1, 2), Complex<double>(-8.55, -8.32));
    ASSERT_COMPLEX_DOUBLE_EQ(pow(z2, -5), Complex<double>(-4.70251762084e-05, -9.79536631956e-06));
}

// TEST pow(complex, real_type)
TEST(Pow, Complex_to_real)
{
    ASSERT_COMPLEX_DOUBLE_EQ(pow(z1, -2.51), Complex<double>(-0.0439236260895, 0.0074249581562));
    ASSERT_COMPLEX_DOUBLE_EQ(pow(z2, 3.1415), Complex<double>(-491.242404853, -160.345650462));
}

// TEST pow(real_type, complex)
TEST(Pow, Real_to_complex)
{
    ASSERT_COMPLEX_DOUBLE_EQ(pow(value, z1), Complex<double>(11.0825338984, -15.2095401654));
    ASSERT_COMPLEX_DOUBLE_EQ(pow(value, z2), Complex<double>(5.42440728125e-08, 4.36212648917e-08));
}

// TEST pow(complex, complex)
TEST(Pow, Complex_to_complex)
{
    ASSERT_COMPLEX_DOUBLE_EQ(pow(z1, z2), Complex<double>(-0.000154301031075, 6.57928589888e-05));
    ASSERT_COMPLEX_DOUBLE_EQ(pow(z2, z1), Complex<double>(-186635.911010941, -195414.098065807));
}

// TEST Trigonometric functions
TEST(Trigonometric, Sin)
{
    ASSERT_COMPLEX_DOUBLE_EQ(sin(z1), Complex<double>(11.8388985179, -3.27575962455));
    ASSERT_COMPLEX_DOUBLE_EQ(sin(z2), Complex<double>(-0.888994153369, 0.160201279969));
}

TEST(Trigonometric, Sinh)
{
    ASSERT_COMPLEX_DOUBLE_EQ(sinh(z1), Complex<double>(-1.69548631445, 0.115050429965));
    ASSERT_COMPLEX_DOUBLE_EQ(sinh(z2), Complex<double>(-707.091945152, 218.729370078));
}

TEST(Trigonometric, aSinh)
{
    ASSERT_COMPLEX_DOUBLE_EQ(asinh(z1), Complex<double>(1.917661428, -1.16961813));
    ASSERT_COMPLEX_DOUBLE_EQ(asinh(z2), Complex<double>(-2.686500788, 0.04069378972));
}

TEST(Trigonometric, aSin)
{
    ASSERT_COMPLEX_DOUBLE_EQ(asin(z1), Complex<double>(0.371904254, -1.947656631));
    ASSERT_COMPLEX_DOUBLE_EQ(asin(z2), Complex<double>(-1.52933383, 2.677164793));
}

TEST(Trigonometric, Cos)
{
    ASSERT_COMPLEX_DOUBLE_EQ(cos(z1), Complex<double>(3.28666346637, 11.7996217626));
    ASSERT_COMPLEX_DOUBLE_EQ(cos(z2), Complex<double>(0.549929090336, 0.258975209272));
}

TEST(Trigonometric, Cosh)
{
    ASSERT_COMPLEX_DOUBLE_EQ(cosh(z1), Complex<double>(-1.96755337967, 0.0991416199894));
    ASSERT_COMPLEX_DOUBLE_EQ(cosh(z2), Complex<double>(707.092590519, -218.729170443));
}

TEST(Trigonometric, aCosh)
{
    ASSERT_COMPLEX_DOUBLE_EQ(acosh(z1), Complex<double>(1.947656631, -1.198892117));
    ASSERT_COMPLEX_DOUBLE_EQ(acosh(z2), Complex<double>(2.677164793, 3.100130081));
}

TEST(Trigonometric, aCos)
{
    ASSERT_COMPLEX_DOUBLE_EQ(acos(z1), Complex<double>(1.198892117, 1.947656631));
    ASSERT_COMPLEX_DOUBLE_EQ(acos(z2), Complex<double>(3.100130081, -2.677164793));
}

TEST(Trigonometric, Tan)
{
    ASSERT_COMPLEX_DOUBLE_EQ(tan(z1), Complex<double>(0.00171795731576, -1.00285012591));
    ASSERT_COMPLEX_DOUBLE_EQ(tan(z2), Complex<double>(-1.21084572844, 0.861529812139));
}

TEST(Trigonometric, Tanh)
{
    ASSERT_COMPLEX_DOUBLE_EQ(tanh(z1), Complex<double>(0.8624796867, -0.01501500234));
    ASSERT_COMPLEX_DOUBLE_EQ(tanh(z2), Complex<double>(-0.9999992847, 5.153517009e-07));
}

TEST(Trigonometric, aTanh)
{
    ASSERT_COMPLEX_DOUBLE_EQ(atanh(z1), Complex<double>(0.1019303277, -1.305935025));
    ASSERT_COMPLEX_DOUBLE_EQ(atanh(z2), Complex<double>(-0.137613073, 1.565069199));
}

TEST(Trigonometric, aTan)
{
    ASSERT_COMPLEX_DOUBLE_EQ(atan(z1), Complex<double>(1.454027772, -0.2713128328));
    ASSERT_COMPLEX_DOUBLE_EQ(atan(z2), Complex<double>(-1.434879899, 0.005516957957));
}

TEST(Trigonometric, Cot)
{
    ASSERT_COMPLEX_DOUBLE_EQ(cot(z1), Complex<double>(0.00170820122099, 0.997155047942));
    ASSERT_COMPLEX_DOUBLE_EQ(cot(z2), Complex<double>(-0.548295665932, -0.390118287552));
}

TEST(Trigonometric, Coth)
{
    ASSERT_COMPLEX_DOUBLE_EQ(coth(z1), Complex<double>(1.159096241, 0.02017883584));
    ASSERT_COMPLEX_DOUBLE_EQ(coth(z2), Complex<double>(-1.000000715, -4.923758752e-07));
}

TEST(Trigonometric, aCoth)
{
    ASSERT_COMPLEX_DOUBLE_EQ(acoth(z1), Complex<double>(0.1019303277, 0.2648612559));
    ASSERT_COMPLEX_DOUBLE_EQ(acoth(z2), Complex<double>(-0.137613073, -0.005727126263));
}

TEST(Trigonometric, aCot)
{
    ASSERT_COMPLEX_DOUBLE_EQ(acot(z1), Complex<double>(0.1167684942, 0.2713128328));
    ASSERT_COMPLEX_DOUBLE_EQ(acot(z2), Complex<double>(-0.135916397, -0.005516957492));
}

TEST(Trigonometric, Sec)
{
    ASSERT_COMPLEX_DOUBLE_EQ(sec(z1), Complex<double>(0.021906236744, -0.0786467219612));
    ASSERT_COMPLEX_DOUBLE_EQ(sec(z2), Complex<double>(1.48834568491, -0.700898792213));
}

TEST(Trigonometric, Sech)
{
    ASSERT_COMPLEX_DOUBLE_EQ(sech(z1), Complex<double>(-0.5069583058, -0.02554477379));
    ASSERT_COMPLEX_DOUBLE_EQ(sech(z2), Complex<double>(0.001290732995, 0.0003992701822));
}

TEST(Trigonometric, aSech)
{
    ASSERT_COMPLEX_DOUBLE_EQ(asech(z1), Complex<double>(0.266560357, 1.46539224358));
    ASSERT_COMPLEX_DOUBLE_EQ(asech(z2), Complex<double>(0.00567335253, -1.7079793389));
}

TEST(Trigonometric, aSec)
{
    ASSERT_COMPLEX_DOUBLE_EQ(asec(z1), Complex<double>(1.4653922436, -0.26656035702));
    ASSERT_COMPLEX_DOUBLE_EQ(asec(z2), Complex<double>(1.7079793389, 0.0056733525294));
}

TEST(Trigonometric, Csc)
{
    ASSERT_COMPLEX_DOUBLE_EQ(csc(z1), Complex<double>(0.0784603960681, 0.0217095701242));
    ASSERT_COMPLEX_DOUBLE_EQ(csc(z2), Complex<double>(-1.08948692501, -0.196331099855));
}

TEST(Trigonometric, Csch)
{
    ASSERT_COMPLEX_DOUBLE_EQ(csch(z1), Complex<double>(-0.5870980024, -0.03983867913));
    ASSERT_COMPLEX_DOUBLE_EQ(csch(z2), Complex<double>(-0.00129073381, -0.0003992710845));
}

TEST(Trigonometric, aCsch)
{
    ASSERT_COMPLEX_DOUBLE_EQ(acsch(z1), Complex<double>(0.11281932446, 0.26979442619));
    ASSERT_COMPLEX_DOUBLE_EQ(acsch(z2), Complex<double>(-0.13633472073, -0.0055682820015));
}

TEST(Trigonometric, aCsc)
{
    ASSERT_COMPLEX_DOUBLE_EQ(acsc(z1), Complex<double>(0.1054040832, 0.266560357032));
    ASSERT_COMPLEX_DOUBLE_EQ(acsc(z2), Complex<double>(-0.13718301215, -0.0056733525294));
}
#pragma once
#include <limits>
#include <type_traits>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>

#define NODISCARD [[nodiscard]]

constexpr static double PI(3.1415926535897932384626433832795029L);

template<class T>
class Complex;

template<class T>
class Complex_base
{
	T m_Real = static_cast<T>(0);
	T m_Imag = static_cast<T>(0);
public:
	using value_t = T;

	constexpr Complex_base() = default;
	constexpr Complex_base(const T& Real, const T& Imag) : m_Real(Real), m_Imag(Imag) {}

	NODISCARD constexpr auto real() const { return this->m_Real; }
	NODISCARD constexpr auto imag() const { return this->m_Imag; }

	T real(const T& NewReal) { return (this->m_Real = NewReal); }
	T imag(const T& NewImag) { return (this->m_Imag = NewImag); }

protected:
	template<class U>
	void addition(const Complex<U>& Rhs)
	{
		static_assert(std::is_convertible_v<U, value_t>, "'U' cannot be converted to 'T'");
		this->m_Real += static_cast<value_t>(Rhs.real());
		this->m_Imag += static_cast<value_t>(Rhs.imag());
	}
	void addition(const T& Rhs)
	{
		this->m_Real += Rhs;
	}

	template<class U>
	void subtraction(const Complex<U>& Rhs)
	{
		static_assert(std::is_convertible_v<U, value_t>, "'U' cannot be converted to 'T'");
		this->m_Real -= static_cast<value_t>(Rhs.real());
		this->m_Imag -= static_cast<value_t>(Rhs.imag());
	}
	void subtraction(const T& Rhs)
	{
		this->m_Real -= Rhs;
	}

	template<class U>
	void multiply(const Complex<U>& Rhs)
	{
		static_assert(std::is_convertible_v<U, value_t>, "'U' cannot be converted to 'T'");
		const auto rightReal = static_cast<value_t>(Rhs.real());
		const auto rightImag = static_cast<value_t>(Rhs.imag());

		auto temp = (this->m_Real * rightReal) - (this->m_Imag * rightImag);
		this->m_Imag = (this->m_Real * rightImag) + (this->m_Imag * rightReal);
		this->m_Real = temp;
	}
	void multiply(const T& Rhs)
	{
		this->m_Real *= Rhs;
		this->m_Imag *= Rhs;
	}

	template<class U>
	void division(const Complex<U>& Rhs)
	{
		static_assert(std::is_convertible_v<U, value_t>, "'U' cannot be converted to 'T'");
		auto rightReal = static_cast<value_t>(Rhs.real());
		auto rightImag = static_cast<value_t>(Rhs.imag());

		const value_t div = (rightReal * rightReal) + (rightImag * rightImag);

		auto temp = ((this->m_Real * rightReal) + (this->m_Imag * rightImag)) / div;
		this->m_Imag = ((rightReal * this->m_Imag) - (this->m_Real * rightImag)) / div;
		this->m_Real = temp;
	}
	void division(const T& Rhs)
	{
		this->m_Real /= Rhs;
		this->m_Imag /= Rhs;
	}

};

template<class T>
struct Complex_traits
{
	NODISCARD constexpr static Complex_base<T> sin(const T& Real, const T& Imag)
	{
		return Complex_base<T>(std::sin(Real) * std::cosh(Imag),
			std::cos(Real) * std::sinh(Imag));
	}

	NODISCARD constexpr static Complex_base<T> cos(const T& Real, const T& Imag)
	{
		return Complex_base<T>(std::cos(Real) * std::cosh(Imag),
			-std::sin(Real) * std::sinh(Imag));
	}

	NODISCARD constexpr static Complex_base<T> tan(const T& Real, const T& Imag)
	{
		const auto _Real = static_cast<T>(2)* Real;
		const auto _Imag = static_cast<T>(2)* Imag;

		const auto div = std::cos(_Real) + std::cosh(_Imag);

		return Complex_base<T>(std::sin(_Real) / div, std::sinh(_Imag) / div);
	}

	NODISCARD constexpr static Complex_base<T> cot(const T& Real, const T& Imag)
	{
		const auto _Real = static_cast<T>(2)* Real;
		const auto _Imag = static_cast<T>(2)* Imag;

		const auto div = std::cosh(_Imag) - std::cos(_Real);

		return Complex_base<T>(std::sin(_Real) / div, -std::sinh(_Imag) / div);
	}

	NODISCARD constexpr static Complex_base<T> sec(const T& Real, const T& Imag)
	{
		const auto first = std::cos(Real) * std::cosh(Imag);
		const auto second = std::sin(Real) * std::sinh(Imag);

		const auto div = std::pow(first, 2) + std::pow(second, 2);

		return Complex_base<T>(first / div, second / div);
	}

	NODISCARD constexpr static Complex_base<T> csc(const T& Real, const T& Imag)
	{
		const auto first = std::sin(Real) * std::cosh(Imag);
		const auto second = std::cos(Real) * std::sinh(Imag);

		const auto div = std::pow(first, 2) + std::pow(second, 2);

		return Complex_base<T>(first / div, -(second / div));
	}

};

template<class T>
class Complex final : public Complex_base<T>
{
public:
	using value_t = T;
	using Base = Complex_base<T>;

	constexpr Complex() = default;
	constexpr Complex(const Base& Value) : Base(Value) {}
	constexpr Complex(const value_t& Real, const value_t& Imag = value_t(0)) : Base({ Real, Imag }) {}

public:
	Complex& operator=(const Complex& Rhs)
	{
		this->real(Rhs.real());
		this->imag(Rhs.imag());
		return *this;
	}
	Complex& operator=(const value_t& Rhs)
	{
		this->real(Rhs);
		this->imag(static_cast<T>(0));
		return *this;
	}

	template<class U>
	Complex& operator+=(const Complex<U>& Rhs)
	{
		this->addition(Rhs);
		return *this;
	}
	Complex& operator+=(const Complex& Rhs)
	{
		this->addition(Rhs);
		return *this;
	}
	Complex& operator+=(const T& Rhs)
	{
		this->addition(Rhs);
		return *this;
	}

	template<class U>
	Complex& operator-=(const Complex<U>& Rhs)
	{
		this->subtraction(Rhs);
		return *this;
	}
	Complex& operator-=(const Complex& Rhs)
	{
		this->subtraction(Rhs);
		return *this;
	}
	Complex& operator-=(const T& Rhs)
	{
		this->subtraction(Rhs);
		return *this;
	}

	template<class U>
	Complex& operator*=(const Complex<U>& Rhs)
	{
		this->multiply(Rhs);
		return *this;
	}
	Complex& operator*=(const Complex& Rhs)
	{
		this->multiply(Rhs);
		return *this;
	}
	Complex& operator*=(const T& Rhs)
	{
		this->multiply(Rhs);
		return *this;
	}

	template<class U>
	Complex& operator/=(const Complex<U>& Rhs)
	{
		this->division(Rhs);
		return *this;
	}
	Complex& operator/=(const Complex& Rhs)
	{
		this->division(Rhs);
		return *this;
	}
	Complex& operator/=(const T& Rhs)
	{
		this->division(Rhs);
		return *this;
	}
};

// Template binary operator+
template<class T>
NODISCARD Complex<T> operator+(const Complex<T>& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Lhs);
	return temp += Rhs;
}
template<class T>
NODISCARD Complex<T> operator+(const Complex<T>& Lhs, const T& Rhs)
{
	Complex<T> temp(Lhs);
	return temp += Rhs;
}
template<class T>
NODISCARD Complex<T> operator+(const T& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Lhs);
	return temp += Rhs;
}

// Template binary operator-
template<class T>
NODISCARD Complex<T> operator-(const Complex<T>& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Lhs);
	return temp -= Rhs;
}
template<class T>
NODISCARD Complex<T> operator-(const Complex<T>& Lhs, const T& Rhs)
{
	Complex<T> temp(Lhs);
	return temp -= Rhs;
}
template<class T>
NODISCARD Complex<T> operator-(const T& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Lhs);
	return temp -= Rhs;
}

// Template binary operator*
template<class T>
NODISCARD Complex<T> operator*(const Complex<T>& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Lhs);
	return temp *= Rhs;
}
template<class T>
NODISCARD Complex<T> operator*(const Complex<T>& Lhs, const T& Rhs)
{
	Complex<T> temp(Lhs);
	return temp *= Rhs;
}
template<class T>
NODISCARD Complex<T> operator*(const T& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Rhs);
	return temp *= Lhs;
}

// Template binary operator/
template<class T>
NODISCARD Complex<T> operator/(const Complex<T>& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Lhs);
	return temp /= Rhs;
}
template<class T>
NODISCARD Complex<T> operator/(const Complex<T>& Lhs, const T& Rhs)
{
	Complex<T> temp(Lhs);
	return temp /= Rhs;
}
template<class T>
NODISCARD Complex<T> operator/(const T& Lhs, const Complex<T>& Rhs)
{
	Complex<T> temp(Rhs);
	return temp /= Lhs;
}

// Template unary operator+
template<class T>
NODISCARD Complex<T> operator+(const Complex<T>& Lhs)
{
	return Lhs;
}

// Template unary operator-
template<class T>
NODISCARD Complex<T> operator-(const Complex<T>& Lhs)
{
	return Complex<T>(-Lhs.real(), -Lhs.imag());
}

// Template binary operator==
template<class T>
NODISCARD constexpr bool operator==(const Complex<T>& Lhs, const Complex<T>& Rhs)
{
	return (Lhs.real() == Rhs.real() && Lhs.imag() == Rhs.imag());
}

template<class T>
NODISCARD constexpr bool operator==(const Complex<T>& Lhs, const T& Rhs)
{
	return (Lhs.real() == Rhs && Lhs.imag() == static_cast<T>(0));
}
template<class T>
NODISCARD constexpr bool operator==(const T& Lhs, const Complex<T>& Rhs)
{
	return (Lhs == Rhs.real() && static_cast<T>(0) == Rhs.imag());
}

// Template binary operator!=
template<class T>
NODISCARD constexpr bool operator!=(const Complex<T>& Lhs, const Complex<T>& Rhs)
{
	return !(Lhs == Rhs);
}
template<class T>
NODISCARD constexpr bool operator!=(const Complex<T>& Lhs, const T& Rhs)
{
	return !(Lhs == Rhs);
}
template<class T>
NODISCARD constexpr bool operator!=(const T& Lhs, const Complex<T>& Rhs)
{
	return !(Lhs == Rhs);
}

// Magnitude function
template<class T>
NODISCARD constexpr inline T abs(const Complex<T>& Lhs)
{
	return std::hypot(Lhs.real(), Lhs.imag());
}

// Phase function
template<class T>
NODISCARD constexpr inline T arg(const Complex<T>& Lhs)
{
	return std::atan2(Lhs.imag(), Lhs.real());
}

// Norm function
template<class T>
NODISCARD constexpr inline T norm(const Complex<T>& Lhs)
{
	return Lhs.real() * Lhs.real() + Lhs.imag() * Lhs.imag();
}

// Conjugation function
template<class T>
NODISCARD constexpr inline Complex<T> conjugation(const Complex<T>& Lhs)
{
	return Complex<T>(Lhs.real(), -Lhs.imag());
}

// Get real part
template<class T>
NODISCARD constexpr inline auto real(const Complex<T>& Lhs)
{
	return Lhs.real();
}

// Get imag part
template<class T>
NODISCARD constexpr inline auto imag(const Complex<T>& Lhs)
{
	return Lhs.imag();
}

// Exp function
template<class T>
NODISCARD constexpr inline Complex<T> exp(const Complex<T>& Lhs)
{
	const T Imag = Lhs.imag();
	const T _Exp = std::exp(Lhs.real());
	return Complex<T>(_Exp * std::cos(Imag), _Exp * std::sin(Imag));
}

// Log function
template<class T>
NODISCARD constexpr inline Complex<T> log(const Complex<T>& Lhs)
{
	return Complex<T>(std::log(abs(Lhs)), arg(Lhs));
}

// Log base 10 function
template<class T>
NODISCARD constexpr inline Complex<T> log10(const Complex<T>& Lhs)
{
	constexpr auto log10_exp(static_cast<T>(0.43429448190325182765112891891660508L));
	return log(Lhs) * log10_exp;
}

// Pow function
template<class T>
NODISCARD constexpr Complex<T> pow(const Complex<T>& Lhs, const int& Rhs)
{
	if (Lhs.imag() == static_cast<T>(0))
		return Complex<T>(std::pow(Lhs.real(), Rhs), 0);
	else
	{
		Complex<T> Temp(Lhs);
		const size_t Count = std::abs(Rhs);
		for (size_t i = 1; i < Count; ++i)
			Temp *= Lhs;

		return Rhs < 0 ? (Complex<T>(1) / Temp) : Temp;
	}
}

template<class T, class U>
NODISCARD constexpr Complex<T> pow(const Complex<T>& Lhs, const U& Rhs)
{
	static_assert(std::is_arithmetic_v<U>, "Template parameter 'U' must be a arithmetic type.");

	const T Module = std::pow(abs(Lhs), static_cast<T>(Rhs));
	const T Arg = arg(Lhs) * static_cast<T>(Rhs);
	return Complex<T>(Module * std::cos(Arg), Module * std::sin(Arg));
}

template<class T, class U>
NODISCARD constexpr Complex<U> pow(const T Lhs, const Complex<U>& Rhs)
{
	static_assert(std::is_arithmetic_v<T>, "Template parameter 'T' must be a arithmetic type.");
	// e^(Rhs * ln(Lhs))
	const auto Left = static_cast<U>(std::log(Lhs));
	return exp(Rhs * Left);
}

template<class T, class U>
NODISCARD constexpr Complex<T> pow(const Complex<T>& Lhs, const Complex<U>& Rhs)
{
	return exp(Rhs * log(Lhs));
}

// Nth root function
template<class T>
NODISCARD constexpr Complex<T> nth_root(const Complex<T>& Lhs,int Rhs)
{
	if (Rhs == 0)
		throw std::invalid_argument("Invalid root argument(root < 2).");

	const T Module = std::pow(abs(Lhs), static_cast<T>(1) / Rhs);
	T Arg = arg(Lhs);
	
	if (Lhs.imag() < 0)
		Arg += static_cast<T>(PI * 2);
	Arg /= Rhs;

	return Complex<T>(Module * std::cos(Arg), Module * std::sin(Arg));
}

// Sqrt function
template<class T>
NODISCARD constexpr inline Complex<T> sqrt(const Complex<T>& Lhs)
{
	const T Module = std::pow(abs(Lhs), static_cast<T>(0.5));
	const T Arg = arg(Lhs) / static_cast<T>(2.0);
	return Complex<T>(Module * std::cos(Arg), Module * std::sin(Arg));
}

// Sine function
template<class T>
NODISCARD constexpr Complex<T> sin(const Complex<T>& Lhs)
{
	return Complex<T>(Complex_traits<T>::sin(Lhs.real(), Lhs.imag()));
}

// Hyperbolic sine funtion
template<class T>
NODISCARD constexpr Complex<T> sinh(const Complex<T>& Lhs)
{
	const auto Temp = Complex_traits<T>::sin(Lhs.imag(), Lhs.real());
	return Complex<T>(Temp.imag(), Temp.real());
}

// Inverse hyperbolic sine funtion
template<class T>
NODISCARD constexpr Complex<T> asinh(const Complex<T>& Lhs)
{
	auto Temp = pow(Lhs, 2) + static_cast<T>(1);
	return log(Lhs + sqrt(Temp));
}

// Inverse sine function
template<class T>
NODISCARD constexpr Complex<T> asin(const Complex<T>& Lhs)
{
	Complex<T> Temp = asinh(Complex<T>(-Lhs.imag(), Lhs.real()));
	return Complex<T>(Temp.imag(), -Temp.real());
}

// Cosine function
template<class T>
NODISCARD constexpr Complex<T> cos(const Complex<T>& Lhs)
{
	return Complex<T>(Complex_traits<T>::cos(Lhs.real(), Lhs.imag()));
}

// Hyperbolic cosine funtion
template<class T>
NODISCARD constexpr Complex<T> cosh(const Complex<T>& Lhs)
{
	const auto Temp = Complex_traits<T>::cos(Lhs.imag(), Lhs.real());
	return Complex<T>(Temp.real(), -Temp.imag());
}

// Inverse hyperbolic cosine funtion
template<class T>
NODISCARD constexpr Complex<T> acosh(const Complex<T>& Lhs)
{
	const auto First = sqrt(Lhs - static_cast<T>(1));
	const auto Second = sqrt(Lhs + static_cast<T>(1));
	return log(Lhs + First * Second);
}

// Inverse cosine function
template<class T>
NODISCARD constexpr Complex<T> acos(const Complex<T>& Lhs)
{
	constexpr auto PIby2 = static_cast<T>(PI / 2);
	Complex<T> Temp = static_cast<T>(PIby2) - asin(Lhs);
	return Complex<T>(Temp.real(), Temp.imag());
}

// Tangent function
template<class T>
NODISCARD constexpr Complex<T> tan(const Complex<T>& Lhs)
{
	return Complex<T>(Complex_traits<T>::tan(Lhs.real(), Lhs.imag()));
}

// Hyperbolic tangent funtion
template<class T>
NODISCARD constexpr Complex<T> tanh(const Complex<T>& Lhs)
{
	const auto Temp = Complex_traits<T>::tan(Lhs.imag(), Lhs.real());
	return Complex<T>(Temp.imag(), Temp.real());
}

// Inverse hyperbolic tangent funtion
template<class T>
NODISCARD constexpr Complex<T> atanh(const Complex<T>& Lhs)
{
	const auto Temp = static_cast<T>(0.5)* log((static_cast<T>(1) + Lhs) / (static_cast<T>(1) - Lhs));
	return Complex<T>(Temp.real(), Temp.imag());
}

// Inverse tangent function
template<class T>
NODISCARD constexpr Complex<T> atan(const Complex<T>& Lhs)
{
	auto Temp = atanh(Complex<T>(-Lhs.imag(), Lhs.real()));
	return Complex<T>(Temp.imag(), -Temp.real());
}

// Cotangent function
template<class T>
NODISCARD constexpr Complex<T> cot(const Complex<T>& Lhs)
{
	return Complex<T>(Complex_traits<T>::cot(Lhs.real(), Lhs.imag()));
}

// Hyperbolic cotangent funtion
template<class T>
NODISCARD constexpr Complex<T> coth(const Complex<T>& Lhs)
{
	const auto Temp = Complex_traits<T>::cot(Lhs.imag(), Lhs.real());
	return Complex<T>(-Temp.imag(), -Temp.real());
}

// Inverse hyperbolic cotangent  funtion
template<class T>
NODISCARD constexpr Complex<T> acoth(const Complex<T>& Lhs)
{
	const auto Temp = static_cast<T>(0.5)* log((Lhs + static_cast<T>(1)) / (Lhs - static_cast<T>(1)));
	constexpr auto PIby2 = static_cast<T>(PI / 2);

	return Complex<T>(Temp.real(), Temp.imag());
}

// Inverse cotangent function
template<class T>
NODISCARD constexpr Complex<T> acot(const Complex<T>& Lhs)
{
	auto Temp = acoth(Complex<T>(-Lhs.imag(), Lhs.real()));
	return Complex<T>(-Temp.imag(), Temp.real());
}

// Secant function
template<class T>
NODISCARD constexpr Complex<T> sec(const Complex<T>& Lhs)
{
	return Complex<T>(Complex_traits<T>::sec(Lhs.real(), Lhs.imag()));
}

// Hyperbolic secant function
template<class T>
NODISCARD constexpr Complex<T> sech(const Complex<T>& Lhs)
{
	const auto Temp = Complex_traits<T>::sec(Lhs.imag(), Lhs.real());
	return Complex<T>(Temp.real(), -Temp.imag());
}

// Inverse hyperbolic secant function
template<class T>
NODISCARD constexpr Complex<T> asech(const Complex<T>& Lhs)
{
	const auto invers = Complex<T>(1) / Lhs;
	const auto first = sqrt(invers - static_cast<T>(1));
	const auto second = sqrt(invers + static_cast<T>(1));

	return log(first * second + invers);
}

// Inverse secant function
template<class T>
NODISCARD constexpr Complex<T> asec(const Complex<T>& Lhs)
{
	auto Temp = asech(Lhs);
	const auto Imag = Lhs.imag() < 0 ? -Temp.real() : Temp.real();
	return Complex<T>(std::abs(Temp.imag()), Imag);
}

// Cosecant function
template<class T>
NODISCARD constexpr Complex<T> csc(const Complex<T>& Lhs)
{
	return Complex<T>(Complex_traits<T>::csc(Lhs.real(), Lhs.imag()));
}

// Hyperbolic cosecant function
template<class T>
NODISCARD constexpr Complex<T> csch(const Complex<T>& Lhs)
{
	const auto Temp = Complex_traits<T>::csc(Lhs.imag(), Lhs.real());
	return Complex<T>(-Temp.imag(), -Temp.real());
}

// Inverse hyperbolic cosecant function
template<class T>
NODISCARD constexpr Complex<T> acsch(const Complex<T>& Lhs)
{
	const auto invers = Complex<T>(1) / Lhs;
	return log(sqrt(pow(invers, 2) + static_cast<T>(1)) + invers);
}

// Inverse cosecant function
template<class T>
NODISCARD constexpr Complex<T> acsc(const Complex<T>& Lhs)
{
	auto Temp = acsch(Complex<T>(Lhs.imag(), Lhs.real()));
	return Complex<T>(-Temp.imag(), -Temp.real());
}

template<class T>
std::ostream& operator<<(std::ostream& out, const Complex<T>& Lhs)
{
	auto Imag = Lhs.imag();
	out << Lhs.real();
	if (Imag < 0)
		out << " - " << static_cast<T>(-1)* Imag;
	else
		out << " + " << Imag;
	out << 'i';
	return out;
}

template<class T>
NODISCARD constexpr inline Complex<T> make_complex_real(T&& Value)
{
	return Complex<T>(std::forward<T>(Value), static_cast<T>(0));
}

template<class T>
NODISCARD constexpr inline Complex<T> make_complex_imag(T&& Value)
{
	return Complex<T>(static_cast<T>(0), std::forward<T>(Value));
}

template<class T>
NODISCARD constexpr inline Complex<T> make_complex(T&& Real, T&& Imag)
{
	return Complex<T>(std::forward<T>(Real), std::forward<T>(Imag));
}
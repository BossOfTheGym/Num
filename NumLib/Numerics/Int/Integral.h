#pragma once

namespace Num
{
	namespace Int
	{
		//lazy version
		template<class Value, class Count, class Function>
		Value intSimpson_1_3(Function&& func, Value a, Value b, Count n)
		{
			const Value h = static_cast<Value>((b - a) / n);

			Value integral = static_cast<Value>(0);

			for (Count i = static_cast<Count>(0); i < n; i++)
			{
				const Value val0 = func(a + h * (i      ));
				const Value val1 = func(a + h * (i + 0.5));
				const Value val2 = func(a + h * (i + 1.0));

				integral += val0 + 4 * val1 + val2;
			}

			return integral * h / 6;
		}


		//lazy version
		template<class Value, class Count, class Function>
		Value intSimpson_3_8(Function&& func, Value a, Value b, Count n)
		{
			const Value h = static_cast<Value>((b - a) / n);

			Value integral = static_cast<Value>(0);

			for (Count i = static_cast<Count>(0); i < n; i++)
			{
				const Value val0 = func(a + h * (i          ));
				const Value val1 = func(a + h * (i + 1.0 / 3));
				const Value val2 = func(a + h * (i + 2.0 / 3));
				const Value val3 = func(a + h * (i + 1.0    ));

				integral += val0 + 3 * val1 + 3 * val2 + val3;
			}

			return integral * h / 8;
		}


		template<class Value, class Count, class Function>
		Value intLeftRectangles(Function& func, Value a, Value b, Count n)
		{
			const Value h = static_cast<Value>((b - a) / n);

			Value integral = static_cast<Value>(0);

			for (Count i = static_cast<Count>(0); i < n; i++)
			{
				const Value val0 = func(a + h * i);

				integral += val0;
			}

			return integral * h;
		}


		template<class Value, class Count, class Function>
		Value intRightRectangles(Function&& func, Value a, Value b, Count n)
		{
			const Value h = static_cast<Value>((b - a) / n);

			Value integral = static_cast<Value>(0);

			for (Count i = static_cast<Count>(0); i < n; i++)
			{
				const Value val0 = func(a + h * (i + 1));

				integral += val0;
			}

			return integral * h;
		}


		template<class Value, class Count, class Function>
		Value intMidpointRectangles(Function&& func, Value a, Value b, Count n)
		{
			const Value h = static_cast<Value>((b - a) / n);

			Value integral = static_cast<Value>(0);

			for (Count i = static_cast<Count>(0); i < n; i++)
			{
				const Value val0 = func(a + h * (i + 1.0 / 2));

				integral += val0;
			}

			return integral * h;
		}


		template<class Value, class Count, class Function>
		Value intTrapezoidalRule(Function&& func, Value a, Value b, Count n)
		{
			const Value h = static_cast<Value>((b - a) / n);

			Value integral = static_cast<Value>(0);

			for (Count i = static_cast<Count>(1); i < n; i++)
			{
				const double val0 = func(a + h * i);

				integral += val0;
			}

			return (integral + (func(a) + func(b)) / 2) * h;
		}


		template<class Value, class Count, class Function>
		Value intGaussRule_2(Function&& func, Value a, Value b, Count n)
		{
			//Gauss args
			const Value X0 = static_cast<Value>(-1.0 / sqrt(3));
			const Value X1 = static_cast<Value>(+1.0 / sqrt(3));

			const Value h = static_cast<Value>((b - a) / n);

			Value integral = static_cast<Value>(0);

			const Value A = static_cast<Value>(h / 2);
			for (Count i = static_cast<int>(0); i < n; i++)
			{
				const double B = a + h * (i + 0.5);

				integral += func(A * X0 + B);
				integral += func(A * X1 + B);
			}

			return integral * A;
		}
	}
}
package main

import (
	"github.com/borgishmorg/go3dplot"
	"math"
)

const (
	stepT = 101
	stepX = 101
	T     = 1
	l     = 1
	tau   = T / float64(stepT-1)
	h     = l / float64(stepX-1)
	alpha = 0.5
)

func derivative(order float64, f []float64) float64 {
	var sum float64
	for i := 0; i < len(f)-1; i++ {
		sum += (math.Pow(float64(stepX*(len(f)-i-1)), 1-order) - math.Pow(float64(stepX*(len(f)-2-i)), 1-order)) * (f[i+1] - f[i])
	}
	sum /= math.Gamma(2-order) * stepX
	return sum
}

func phi(x float64) float64 {
	return (1 - x) * x
}

func eta(x float64) float64 {
	return (1 - x) * x
}

func fTruth() [stepT][stepX]float64 {
	f := [stepT][stepX]float64{}
	var x, t float64
	for i := 2; i < stepT-1; i++ {
		for j := 1; j < stepX; j++ {
			t = float64(i) * tau
			x = float64(j) * h
			f[i][j] = math.Exp(t) * ((1-x)*x - (1-4*x)/math.Gamma(0.5)/math.Sqrt(x)) // p0 = 1000 * t^2 * (1 - t)
		}
	}
	return f
}

func Gauss(A [stepX - 2][stepX - 2]float64, b [stepX - 2]float64) [stepX - 2]float64 {
	var answer, buff1 [stepX - 2]float64
	var sum, buff2 float64
	for k := 0; k < len(A); k++ {
		max := -100000.0
		num := -1
		// ищем ведущий элемент
		for i := k; i < len(A); i++ {
			if math.Abs(A[i][k]) > max {
				max = math.Abs(A[i][k])
				num = i
			}
		}
		// меняем местами строки
		if num != k {
			buff2 = b[num]
			buff1 = A[num]
			b[num] = b[k]
			A[num] = A[k]
			A[k] = buff1
			b[k] = buff2
		}
		// делим строку на ведущий элемент
		for i := k + 1; i < len(A); i++ {
			A[k][i] /= A[k][k]
		}
		b[k] /= A[k][k]
		A[k][k] = 1
		// вычитаем строку из остальных
		for i := k + 1; i < len(A); i++ {
			for j := k + 1; j < len(A); j++ {
				A[i][j] -= A[k][j] * A[i][k]
			}
			b[i] -= b[k] * A[i][k]
			A[i][k] = 0
		}
	}
	// обратный ход
	for i := len(A) - 1; i >= 0; i-- {
		sum = b[i]
		for j := len(A) - 1; j > i; j-- {
			sum -= A[i][j] * answer[i+1]
		}
		answer[i] = sum
	}
	return answer
}

func implicitStraight(f [stepT][stepX]float64) {
	U := make([][]float64, 2)
	x := make([]float64, 0)
	t := make([]float64, 0)
	var buff float64
	for j := 0; j < stepX; j++ {
		U[0] = append(U[0], phi(h*float64(j)))
		U[1] = append(U[1], U[0][j]+tau*eta(h*float64(j)))
	}
	coefficient := 1 / (math.Gamma(1-alpha) * h * h)
	for n := 2; n < stepT; n++ {
		U = append(U, make([]float64, 0))
		U[n] = append(U[n], 0)
		//создание и заполнение матрицы
		var A [stepX - 2][stepX - 2]float64
		var columnFreeMembers [stepX - 2]float64
		for i := 1; i < stepX-1; i++ {
			A[i-1][i-1] += 1 / (tau * tau)
			columnFreeMembers[i-1] += (2*U[n-1][i]-U[n-2][i])/(tau*tau) + f[n][i]
			for k := 1; k <= i+1; k++ {
				buff = coefficient * (math.Pow(h*float64(i+1-k), 1-alpha)*h*(alpha-3+float64(k-i)) + math.Pow(h*float64(i+2-k), 2-alpha)) / (h * (alpha*alpha - 3*alpha + 2))
				if i < stepX-2 || k != i+1 {
					A[i-1][k-1] -= buff
				}
				if k-1 > 0 {
					A[i-1][k-2] += buff - coefficient*(math.Pow(h*float64(i+1-k), 1-alpha)-math.Pow(h*float64(i+2-k), 1-alpha))/(alpha-1)
				}
			}
			for k := 1; k <= i; k++ {
				buff = 2 * coefficient * (math.Pow(h*float64(i-k), 1-alpha)*h*(alpha-2+float64(k-i)) + math.Pow(h*float64(i+1-k), 2-alpha)) / (h * (alpha*alpha - 3*alpha + 2))
				A[i-1][k-1] += buff
				if k-1 > 0 {
					A[i-1][k-2] -= buff - 2*coefficient*(math.Pow(h*float64(i-k), 1-alpha)-math.Pow(h*float64(i+1-k), 1-alpha))/(alpha-1)
				}
			}
			for k := 1; k <= i-1; k++ {
				buff = coefficient * (math.Pow(h*float64(i-1-k), 1-alpha)*h*(alpha-1+float64(k-i)) + math.Pow(h*float64(i-k), 2-alpha)) / (h * (alpha*alpha - 3*alpha + 2))
				A[i-1][k-1] -= buff
				if k-1 > 0 {
					A[i-1][k-2] += buff - coefficient*(math.Pow(h*float64(i-1-k), 1-alpha)-math.Pow(h*float64(i-k), 1-alpha))/(alpha-1)
				}
			}
		}
		//решаем систему
		answer := Gauss(A, columnFreeMembers)
		//заполняем слой
		for j := 0; j < len(answer); j++ {
			U[n] = append(U[n], answer[j])
		}
		U[n] = append(U[n], 0)
	}
	for i := 0; i < stepX; i++ {
		x = append(x, float64(i)*h)
	}
	for i := 0; i < stepT; i++ {
		t = append(t, float64(i)*tau)
	}
	drawer := go3dplot.GetGnuplotDrawer()
	err := drawer.Draw(x, t, U, "example")
	if err != nil {
		return
	}
}

func main() {
	f := fTruth()
	x := make([]float64, 0)
	t := make([]float64, 0)
	U := make([][]float64, 0)
	implicitStraight(f)
	for i := 0; i < stepX; i++ {
		x = append(x, float64(i)*h)
	}
	for i := 0; i < stepT; i++ {
		t = append(t, float64(i)*tau)
	}
	for i := 0; i < stepT; i++ {
		U = append(U, make([]float64, 0))
		for j := 0; j < stepX; j++ {
			U[i] = append(U[i], math.Exp(float64(i)*tau)*((1-float64(j)*h)*float64(j)*h))
		}
	}
	drawer := go3dplot.GetGnuplotDrawer()
	err := drawer.Draw(x, t, U, "example1")
	if err != nil {
		return
	}
}

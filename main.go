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

func phi(x float64) float64 {
	return (1 - x) * x
}

func eta(x float64) float64 {
	return (1 - x) * x
}

func impact(f *[stepT][stepX]float64) {
	var x, t float64
	for i := 2; i < stepT; i++ {
		for j := 1; j < stepX-1; j++ {
			t = float64(i) * tau
			x = float64(j) * h
			f[i][j] = math.Exp(t) * ((1-x)*x - (math.Pow(x, -alpha)*(alpha+2*x-1)/(alpha-1))/math.Gamma(1-alpha))
		}
	}
}

func gauss(matrix *[stepX - 2][stepX - 2]float64, columnFreeMembers *[stepX - 2]float64) [stepX - 2]float64 {
	var answer, buffer1 [stepX - 2]float64
	var buffer2 float64
	for k := 0; k < len(matrix); k++ {
		max := -1e5
		number := -1
		// ищем ведущий элемент
		for i := k; i < len(matrix); i++ {
			if math.Abs(matrix[i][k]) > max {
				max = math.Abs(matrix[i][k])
				number = i
			}
		}
		// меняем местами строки
		if number != k {
			buffer1, buffer2 = matrix[number], columnFreeMembers[number]
			matrix[number], columnFreeMembers[number] = matrix[k], columnFreeMembers[k]
			matrix[k], columnFreeMembers[k] = buffer1, buffer2
		}
		// делим строку на ведущий элемент
		for i := k + 1; i < len(matrix); i++ {
			matrix[k][i] /= matrix[k][k]
		}
		columnFreeMembers[k] /= matrix[k][k]
		matrix[k][k] = 1
		// вычитаем строку из остальных
		for i := k + 1; i < len(matrix); i++ {
			for j := k + 1; j < len(matrix); j++ {
				matrix[i][j] -= matrix[k][j] * matrix[i][k]
			}
			columnFreeMembers[i] -= columnFreeMembers[k] * matrix[i][k]
			matrix[i][k] = 0
		}
	}
	// обратный ход
	for i := len(matrix) - 1; i >= 0; i-- {
		buffer2 = columnFreeMembers[i]
		for j := len(matrix) - 1; j > i; j-- {
			buffer2 -= matrix[i][j] * answer[i+1]
		}
		answer[i] = buffer2
	}
	return answer
}

func implicitStraight(f *[stepT][stepX]float64) [][]float64 {
	U := make([][]float64, 2)
	var buffer1 float64
	for j := 0; j < stepX; j++ {
		U[0] = append(U[0], phi(h*float64(j)))
		U[1] = append(U[1], U[0][j]+tau*eta(h*float64(j)))
	}
	coefficient := 1 / (math.Gamma(1-alpha) * h * h)
	for n := 2; n < stepT; n++ {
		U = append(U, make([]float64, 0))
		U[n] = append(U[n], 0)
		//создание и заполнение матрицы
		var matrix [stepX - 2][stepX - 2]float64
		var columnFreeMembers [stepX - 2]float64
		for i := 1; i < stepX-1; i++ {
			matrix[i-1][i-1] += 1 / (tau * tau)
			columnFreeMembers[i-1] += (2*U[n-1][i]-U[n-2][i])/(tau*tau) + f[n][i]
			for k := 1; k <= i+1; k++ {
				buffer1 = coefficient * (math.Pow(h*float64(i+1-k), 1-alpha)*h*(alpha-3+float64(k-i)) + math.Pow(h*float64(i+2-k), 2-alpha)) / (h * (alpha*alpha - 3*alpha + 2))
				if i < stepX-2 || k != i+1 {
					matrix[i-1][k-1] -= buffer1
				}
				if k-1 > 0 {
					matrix[i-1][k-2] += buffer1 - coefficient*(math.Pow(h*float64(i+1-k), 1-alpha)-math.Pow(h*float64(i+2-k), 1-alpha))/(alpha-1)
				}
			}
			for k := 1; k <= i; k++ {
				buffer1 = 2 * coefficient * (math.Pow(h*float64(i-k), 1-alpha)*h*(alpha-2+float64(k-i)) + math.Pow(h*float64(i+1-k), 2-alpha)) / (h * (alpha*alpha - 3*alpha + 2))
				matrix[i-1][k-1] += buffer1
				if k-1 > 0 {
					matrix[i-1][k-2] -= buffer1 - 2*coefficient*(math.Pow(h*float64(i-k), 1-alpha)-math.Pow(h*float64(i+1-k), 1-alpha))/(alpha-1)
				}
			}
			for k := 1; k <= i-1; k++ {
				buffer1 = coefficient * (math.Pow(h*float64(i-1-k), 1-alpha)*h*(alpha-1+float64(k-i)) + math.Pow(h*float64(i-k), 2-alpha)) / (h * (alpha*alpha - 3*alpha + 2))
				matrix[i-1][k-1] -= buffer1
				if k-1 > 0 {
					matrix[i-1][k-2] += buffer1 - coefficient*(math.Pow(h*float64(i-1-k), 1-alpha)-math.Pow(h*float64(i-k), 1-alpha))/(alpha-1)
				}
			}
		}
		//решаем систему
		answer := gauss(&matrix, &columnFreeMembers)
		//заполняем слой
		for j := 0; j < len(answer); j++ {
			U[n] = append(U[n], answer[j])
		}
		U[n] = append(U[n], 0)
	}
	return U
}

func main() {
	var f [stepT][stepX]float64
	impact(&f)
	x := make([]float64, 0)
	t := make([]float64, 0)
	U := implicitStraight(&f)
	for i := 0; i < stepX; i++ {
		x = append(x, float64(i)*h)
	}
	for i := 0; i < stepT; i++ {
		t = append(t, float64(i)*tau)
	}
	//UErr := make([][]float64, 0)
	//UT := make([][]float64, 0)
	//for i := 0; i < stepT; i++ {
	//	UT = append(UT, make([]float64, 0))
	//	for j := 0; j < stepX; j++ {
	//		UT[i] = append(UT[i], math.Exp(float64(i)*tau)*((1-float64(j)*h)*float64(j)*h))
	//	}
	//}
	//for i := 0; i < stepT; i++ {
	//	UErr = append(UErr, make([]float64, 0))
	//	for j := 0; j < stepX; j++ {
	//		UErr[i] = append(UErr[i], math.Abs(UT[i][j]-U[i][j]))
	//	}
	//}
	// отрисовка графика функции U(x,t)
	drawer := go3dplot.GetGnuplotDrawer()
	err := drawer.Draw(x, t, U, "example")
	if err != nil {
		return
	}
}

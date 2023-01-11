package main

import (
	"fmt"
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
	return math.Exp(x) * (1 - x) * x
}

func eta(x float64) float64 {
	return x * 0
}

func fTruth() [stepT][stepX]float64 {
	f := [stepT][stepX]float64{}
	var x, t float64
	for i := 0; i < stepT; i++ {
		for j := 0; j < stepX; j++ {
			t = float64(i) * tau
			x = float64(j) * h
			f[i][j] = 0 * x * t // p0 = 1000 * t^2 * (1 - t)
		}
	}
	return f
}

func diagonal(A [][]float64) []float64 {
	mtr := MakeRatio(A)
	answer := make([]float64, 0)
	answer = append(answer, mtr[1][len(A)-1])
	for i := 2; i <= len(A); i++ {
		answer = append(answer, mtr[0][len(A)-i]*answer[i-2]+mtr[1][len(A)-i])
	}
	return answer
}

func MakeRatio(A [][]float64) [][]float64 {
	mtr := make([][]float64, 0)
	var v float64
	alpha0 := make([]float64, 0)
	beta := make([]float64, 0)
	alpha0 = append(alpha0, -A[0][2]/A[0][1])
	beta = append(beta, A[0][3]/A[0][1])
	for i := 1; i < len(A); i++ {
		v = -A[i][1] - A[i][0]*alpha0[i-1]
		if i == len(A)-1 {
			alpha0 = append(alpha0, 0)
		} else {
			alpha0 = append(alpha0, A[i][2]/v)
		}
		beta = append(beta, (A[i][0]*beta[i-1]-A[i][3])/v)
	}
	mtr = append(mtr, alpha0)
	mtr = append(mtr, beta)
	return mtr
}

func implicitStraight(f [stepT][stepX]float64) {
	U := make([][]float64, 2)
	x := make([]float64, 0)
	t := make([]float64, 0)
	for j := 0; j < stepX; j++ {
		U[0] = append(U[0], phi(h*float64(j)))
		U[1] = append(U[1], U[0][j]+tau*eta(h*float64(j)))
	}
	buff := 1 / (math.Gamma(1-alpha) * (alpha*alpha - 3*alpha + 2) * math.Pow(h, 1+alpha))
	for i := 2; i < stepT; i++ {
		U = append(U, make([]float64, 0))
		U[i] = append(U[i], 0)
		//создание и заполнение матрицы
		A := make([][]float64, 0)
		for j := 1; j < stepX-1; j++ {
			A = append(A, make([]float64, 4))
			if j > 1 {
				A[j-1][0] = -buff * (alpha - 1)
			}
			A[j-1][1] = 1/(tau*tau) + alpha*buff
			A[j-1][3] = (2*U[i-1][j]-U[i-2][j])/(tau*tau) + f[i][j]
			if j < stepX-2 {
				A[j-1][2] = -buff
			}
		}
		A[0][0], A[len(A)-1][2] = 0, 0
		//решаем систему
		answer := diagonal(A)
		//заполняем слой
		for j := 1; j <= len(answer); j++ {
			U[i] = append(U[i], answer[len(answer)-j])
		}
		U[i] = append(U[i], 0)
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
	fmt.Print(math.Gamma(1 - alpha))
	implicitStraight(f)
	//mtr := [stepX - 2][4]float64{{0, 2, 3, 1}, {1, 2, 3, 1}, {1, 2, 3, 1}, {1, 2, 0, 1}}
	//answer := diagonal(mtr)
	//fmt.Print(answer)
}

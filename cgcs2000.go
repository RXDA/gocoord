package gocoord

import "math"

// CGCS2000toWGS84
/**
  将大地2000转为WGS84
  高斯投影反算为大地平面。
  x，y 高斯平面坐标点
  L0 通过经纬度来获取中央带所在带的角度
  return B纬度 , L经度
*/
func CGCS2000toWGS84(x, y, l0 float64) (float64, float64) {
	//中央子午线经度
	//WGS-84   椭球体参数
	a := 6378137.0               //major semi axis
	efang := 0.0066943799901413  //square of e
	e2fang := 0.0067394967422764 //suqre of e2
	y -= 500000

	// 主曲率计算
	m0 := a * (1 - efang)
	m2 := 3.0 / 2.0 * efang * m0
	m4 := efang * m2 * 5.0 / 4.0
	m6 := efang * m4 * 7.0 / 6.0
	m8 := efang * m6 * 9.0 / 8.0

	// 子午线曲率计算
	a0 := m0 + m2/2.0 + m4*3.0/8.0 + m6*5.0/16.0 + m8*35.0/128.0
	a2 := m2/2.0 + m4/2.0 + m6*15.0/32.0 + m8*7.0/16.0
	a4 := m4/8.0 + m6*3.0/16.0 + m8*7.0/32.0
	a6 := m6/32.0 + m8/16.0
	a8 := m8 / 128.0

	X := x
	FBf := 0.0
	Bf0 := X / a0
	Bf1 := 0.0

	//计算Bf的值，直到满足条件
	for Bf0-Bf1 >= 0.0001 {
		Bf1 = Bf0
		FBf = -a2*math.Sin(2*Bf0)/2 + a4*math.Sin(4*Bf0)/4 - a6*math.Sin(6*Bf0)/6 + a8*math.Sin(8*Bf0)/8
		Bf0 = (X - FBf) / a0
	}
	Bf := Bf0
	//计算公式中参数
	Wf := math.Sqrt(1 - efang*math.Sin(Bf)*math.Sin(Bf))
	Nf := a / Wf
	Mf := a * (1 - efang) / math.Pow(Wf, 3)
	nffang := e2fang * math.Pow(math.Cos(Bf), 2)
	tf := math.Tan(Bf)
	B := Bf - tf*y*y/(2*Mf*Nf) + tf*(5+3*tf*tf+nffang-9*nffang*tf*tf)*math.Pow(y, 4)/(24*Mf*math.Pow(Nf, 3)) -
		tf*(61+90*tf*tf+45*math.Pow(tf, 4))*math.Pow(y, 6)/(720*Mf*math.Pow(Nf, 5))
	l := y/(Nf*math.Cos(Bf)) - (1+2*tf*tf+nffang)*math.Pow(y, 3)/(6*math.Pow(Nf, 3)*math.Cos(Bf)) +
		(5+28*tf*tf+24*math.Pow(tf, 4))*math.Pow(y, 5)/(120*math.Pow(Nf, 5)*math.Cos(Bf))
	L := l + l0
	Bdec := rad2dec(B)
	Ldec := rad2dec(L)
	//转化成为十进制经纬度格式
	return Bdec, Ldec
}

func rad2dec(rad float64) float64 {
	p := 180.0 / math.Pi * 3600
	dms := rad * p
	a0 := math.Floor(dms / 3600.0)
	a1 := math.Floor((dms - a0*3600) / 60.0)
	a2 := float64(int64(math.Floor(dms-a0*3600))) - a1*60
	return a0 + a1/60.0 + a2/3600.0
}

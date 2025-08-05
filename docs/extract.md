### `EXTRACTION` 子例程主要过程总结

`EXTRACTION` 子例程的主要目的是计算植物根系从土壤中提取水分的过程，同时考虑了土壤水分含量、地下水位、根系活性等因素。以下是该子例程的主要步骤及对应的公式总结：

#### 1. 初始化参数

设定一些常量，如叶片水势 `potleaf`、萎蔫点水势 `potwilt` 和田间持水量水势 `potfc`，并初始化一些变量。

#### 2. 确定地下水位位置

通过循环找到地下水位所在的层 `iwtd` 和其上方的层 `kwtd`。

#### 3. 计算根系活动层

找到根系活动的最低层 `kroot`。

#### 4. 计算各层的提取易度函数 `easy`

1. 计算每层的中点位置 `vctr4(k)`：
   - \(vctr4(k) = 0.5\times(slz(k)+slz(k + 1))\)
2. 根据土壤深度选择土壤类型 `nsoil`：
   - 若 \(slz(k)< - 0.30\)，则 \(nsoil = soiltxt(1)\)
   - 否则，\(nsoil = soiltxt(2)\)
3. 计算饱和含水量 `smoisat` 和饱和水势 `psisat`：
   - \(smoisat=\theta_{sat}(nsoil)\times\max(\min(\exp((vctr4(k)+1.5)/fdepth),1.),0.1)\)
   - \(psisat = slpots(nsoil)\times\min(\max(\exp(-(vctr4(k)+1.5)/fdepth),1.),10.)\)
4. 计算土壤水势 `pot`：
   - \$(pot = psisat\times(\frac{smoisat}{smoi(k)})^{slbs(nsoil)}\)$
5. 计算土壤因子 `soilfactor`：
   - 若 \(icefac(k)=0\)，则 \(soilfactor = 1\)
   - 否则，\(soilfactor = 0\)
6. 计算提取易度函数 `easy(k)`：
   - \(easy(k)=\max(-\frac{(potleaf - pot)\times soilfactor}{hveg - vctr4(k)},0)\)

#### 5. 消除小的根系活动

计算活动层的最大提取易度 `maxeasy`，并将小于 `0.001 * maxeasy` 的 `easy` 值设为 0。

#### 6. 计算根系活性 `rootactivity`

1. 计算总提取易度 `toteasy`：
   - \(toteasy=\sum_{k = 1}^{nzg}easy(k)\times dz2(k)\)
2. 若 \(toteasy = 0\)，则 \(rootactivity = 0\)
   - 否则，\(rootactivity=\min(\max(\frac{easy\times dz2}{toteasy},0),1)\)

#### 7. 更新非活动天数 `inactivedays`

若 `easy(k) = 0`，则 \(inactivedays(k)=inactivedays(k)+1\)
否则，\(inactivedays(k)=0\)

#### 8. 计算根区土壤水分和田间持水量

1. 计算每层的最小含水量 `smoimin` 和田间持水量 `smoifc`：
   - \(smoimin = smoisat\times(\frac{psisat}{potwilt})^{\frac{1}{slbs(nsoil)}}\)
   - \(smoifc = smoisat\times(\frac{psisat}{potfc})^{\frac{1}{slbs(nsoil)}}\)
2. 计算每层可提取的最大水量 `maxwat(k)`：
   - \(maxwat(k)=\max((smoi(k)-smoimin)\times dz(k),0)\)
3. 计算根区土壤水分 `rootsmoi` 和根区田间持水量 `rootfc`：
   - \(rootsmoi=\sum_{k = \max(kwtd,1)}^{nzg}\max(rootactivity(k)\times(smoi(k)-smoimin),0)\)
   - \(rootfc=\sum_{k = \max(kwtd,1)}^{nzg}\max(rootactivity(k)\times(smoifc - smoimin),0)\)

#### 9. 计算气孔阻力和蒸散量

1. 计算土壤水分胁迫因子 `fswp`：
   - 若 \(rootsmoi\leq0\)，则 \(fswp = 0\)
   - 若 \(\frac{rootsmoi}{rootfc}\leq1\)，则 \(fswp=\frac{rootsmoi}{rootfc}\)
   - 否则，\(fswp = 1\)
2. 计算冠层气孔阻力 `rs_c`：
   - 若 \(fswp = 0\)，则 \(rs_c = 5000\)
   - 否则，\(rs_c=\min(\frac{rs_c\_factor}{fswp},5000)\)
3. 计算土壤表面阻力 `rs_s`：
   - \(rs_s = 33.5+3.5\times(\frac{\theta_{sat}(nsoil)}{smoi(nzg)})^{2.38}\)
4. 计算冠层和土壤的辐射阻力 `R_c` 和 `R_s`：
   - \(R_c=(delta + gamma)\times ra_c+gamma\times rs_c\)
   - \(R_s = R_s+gamma\times rs_s\)
5. 计算冠层和土壤的蒸散系数 `C_c` 和 `C_s`：
   - \(C_c=\frac{1}{1+\frac{R_a\times R_c}{R_s\times(R_c + R_a)}}\)
   - \(C_s=\frac{1}{1+\frac{R_a\times R_s}{R_c\times(R_s + R_a)}}\)
6. 若 \(lai<0.001\)，则 \(C_c = 0\)
7. 计算冠层和土壤的潜在蒸散量 `pet` 和 `pet_s`：
   - \(pet = C_c\times\frac{petfactor\_c}{delta + gamma\times(1+\frac{rs_c}{ra_a + ra_c})}\)
   - \(pet=\max(\frac{deltat\times pet}{\lambda},0)\)
   - \(pet_s = C_s\times\frac{petfactor\_s}{delta + gamma\times(1+\frac{rs_s}{ra_a + ra_c})}\)
   - \(pet_s=\max(\frac{deltat\times pet_s}{\lambda},0)\)
8. 计算蒸腾水量 `transpwater`：
   - \(transpwater = pet\times10^{-3}\)

#### 10. 提取水分

1. 若 \(toteasy = 0\)，则 \(watdef = transpwater\) 并返回
2. 计算每层提取的水量 `extract`：
   - \(extract=\max(rootactivity(k)\times transpwater,0)\)
3. 若 \(extract\leq maxwat(k)\)，则 \(dsmoi(k)=extract\)
   - 否则，\(dsmoi(k)=maxwat(k)\)，\(watdef = watdef+(extract - maxwat(k))\)

#### 11. 检查结果

确保 \(dsmoi\geq0\)，并检查 \(|watdef - transpwater|\) 是否大于 \(10^9\)，若大于则输出错误信息。

以上就是 `EXTRACTION` 子例程的主要过程及对应的公式总结。

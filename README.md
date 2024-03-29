# APPROXIMATION of solar flares

The previously prepared spectra of the considered solar flare were fitted in the 
XSPEC package with a time resolution of 256 ms and 1024 ms. 
Fitting was performed in Python using the PyXspec, AstroPy, SciPy, MatPlotLib libraries. 
Two phenomenological models were used: a power law and a and a broken power law model. 
The preferred model was selected using the Fisher criterion. 

Исследование тонкой временной структуры солнечных вспышек в жёстком рентгеновском диапазоне по данным эксперимента KONUS-WIND

Основной файл build_models.py _в папке fitting_ - запускается в программе XSPEC

## Описание
Целью данной работы является анализ спектров солнечных вспышек, зарегистрированных с высоким временным разрешением, и аппроксимация спектров. 
Исследование проводится по данным эксперимента Конус-Винд.

Фитирование подготовленных спектров проводилось в программе XSPEC с помощью библиотеки PyXspec. Для каждого спектра строились две феноменологические модели: одностепенная (npow) и двухстепенная (nbknpow). Вычисления и построения графиков выполнены в среде PyCharm, с использованием библиотек AstroPy, SciPy, MatPlotLib.

Основные файлы в папке "fitting"
- _run_fit.py_ - аппроксимация спектров в программе XSPEC
- _build_models.py_ - построение моделей

Входные данные:

_папка data_

Данные для вспышки класса М2.1, произошедшей 26 апреля 2003, 
00:55:01 (03301):
- _KW20150507_T45695_1_bg.pha_ -	файл со спектром солнечной вспышки;
- _KW20150507_T45695_1.rmf_ -	файл с фоном;
- _KW20150507_T45695.arf_ -	файлы с откликом детектора.

![image](https://user-images.githubusercontent.com/62285192/222533877-567b6300-5b33-400a-8e47-824508f0c7fd.png)


Выполнена аппроксимация спектров вспышки класса М2.1, произошедшей 26 апреля 2003, 00:55:01 (03301) в программе XSPEC феноменологическими моделями — одностепенной и двустепенной с изломом.
1.	Одно степенная модель (npow)

![image](https://user-images.githubusercontent.com/62285192/222535482-900b6888-8dbc-4002-b550-dd8026c1a884.png),

где  ![image](https://user-images.githubusercontent.com/62285192/222535509-766bdf4f-7e03-4d21-8c4f-e90314ee926b.png) - показатель степени

2.	Двухстепенная модель с изломом (nbknpow) 

![image](https://user-images.githubusercontent.com/62285192/222535622-9899b911-6be6-4089-b391-4389f64d505c.png),

здесь ![image](https://user-images.githubusercontent.com/62285192/222535536-917fd544-3fcc-4ea2-9529-73787c39162c.png) - значение энергии излома, 
 ![image](https://user-images.githubusercontent.com/62285192/222535571-142acd73-b886-4840-818e-6d17bc03d355.png) - показатель степени до излома, ![image](https://user-images.githubusercontent.com/62285192/222535851-a13c124c-e673-4593-a19b-2195054c8ee8.png) - после излома.



После оценки качества аппроксимации выявлена необходимость проведения исследование спектров в другом расширении. Для этого были определены времена накопления спектров и выявлены спектры с наименьшим разрешением. 
С использованием утилиты SumKonusSpectraRate.exe выполнено суммирование спектров таким образом, чтобы разрешение итоговых спектров было не менее 512 миллисекунд.

Подготовленные спектры профитированы в программе XSPEC с помощью библиотеки Pyxspec. 
В случаях со значительным отклонением от модели были исключены плохие каналы.
По итогу работы для каждого спектра построены графики моделей и по критерию Фишера определена предпочтительная модель.

Одностепенная модель (npow) пример для седьмого спектра:

![image](https://user-images.githubusercontent.com/62285192/222534393-eb21b661-6396-4a52-9be2-1e9255702615.png)


Двухстепенная модель (nbknpow):

![image](https://user-images.githubusercontent.com/62285192/222534529-b8389de1-8c5c-40ba-a09d-b8cbcdb2262b.png)


Получена временная зависимость для основных параметров излучения вспышки М2 2003-04-26 00:55:01 на временном масштабе 1024 мс

![image](https://user-images.githubusercontent.com/62285192/222534819-d5c05583-8a78-4c3d-bef4-616996c4f9ec.png)




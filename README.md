# 说明
基于能态密度泛函的载流子迁移发射模型（Carrier-Migration-and-Emission，CME）

该模型首先需要求解DFT，将DFT求解得到的能态密度函数作为输入；

其中test_funCME_P为主函数，用来研究脉冲产生电子数随光强变化的曲线，坐标皆取对数坐标，test_func为调用函数。

模型优点在于，比起TDSE可以模拟多光子吸收的情况，比起TDDFT计算更为迅速，同时可以计算多发脉冲相互作用。缺点在于只能用作半定量分析。

以下是算法实现流程，用于求解模型

![算法流程图](https://user-images.githubusercontent.com/61829316/221743216-d82d8d2b-906f-4a3e-a6bc-8ccf3a2807d1.png)

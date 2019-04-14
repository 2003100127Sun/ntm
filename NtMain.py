import numpy as np
import pandas as pd
import seaborn as sns
import rpy2.robjects as rob
from rpy2.robjects.packages import importr
stats = importr('stats')
import matplotlib.pyplot as plt
from NewtonDividedDifference import newtonDividedDifference
from StatisticalTest import statisticalTest
from Curve import curve
from ProbabilityDistribution import probabilityDistribution
from Interpolation import interpolation


class ntBoard(object):

    def __init__(self, ):
        pass

    def simulateBatchCurve(self, ):
        a = rob.r(
            '''
            
            # This function will implement a task of drawing curve.
            batch_curve=function(){
                upper=matrix(nrow = 10, ncol = 501)
                mean=matrix(nrow = 10, ncol = 501)
                lower=matrix(nrow = 10, ncol = 501)
                for(i in 1:10){
                    for(j in 1:501){
                        step=0.02*(j-1)
                        upper[i,j]=1/(1+exp(-(2+0.3*i)*(step-2)))
                        mean[i,j]=1/(1+exp(-(2.3+0.3*i)*(step-5)))
                        lower[i,j]=1/(1+exp(-(2.3+0.3*i)*(step-8)))
                    }
                }
                x=seq(0,10,by=0.02)
                # windows()
                plot(x,upper[1,],type = "l", col="slateblue3",xlab = "Time", ylab = "Growth")
                for(k in 1:10){
                    lines(x,upper[k,],type = "l", col="slateblue3")
                    lines(x,mean[k,],type = "l", col="gray24")
                    lines(x,lower[k,],type = "l", col="indianred1")
                    par(new=TRUE)
                }

                # sample points from matrix
                sample_points_raw=c(1,2,3,4,5,6,7,8,9,10)
                sample_points_col=c(51,101,151,201,251,301,351,401,451,501)
                upper_points=matrix(nrow = 10, ncol = 10)
                mean_points=matrix(nrow = 10, ncol = 10)
                lower_points=matrix(nrow = 10, ncol = 10)
                for(p in sample_points_raw){
                    v=1
                    for(q in sample_points_col){
                        upper_points[p,v]=upper[p,q]
                        mean_points[p,v]=mean[p,q]
                        lower_points[p,v]=lower[p,q]
                        v=v+1
                    }
                }
                setwd("C:/Users/Chris_Sun/Desktop/")
                write.csv(upper_points,"upper_points.csv",row.names=FALSE)
                write.csv(mean_points,"mean_points.csv",row.names=FALSE)
                write.csv(lower_points,"lower_points.csv",row.names=FALSE)
                return (upper[1,451])
            }
            
            batch_curve()
            
            '''
        )
        return a

    def mouseCaseStudy(self, ):
        """
        Data pooled from SM/J and 10 LG/J 535 F2 animals by James M. Cheverud.

        .. @Cite:
           ------
           James M. Cheverud. Quantitative Trait Loci for Murine Growth. Genetics, 1996.

        :return:
        """
        weeks = rob.IntVector(np.arange(10) + 1)
        weights = rob.FloatVector([4.83, 8.26, 12.56, 19.21, 24.80, 27.37, 29.72, 31.44, 33.71, 35.51])
        dd = newtonDividedDifference().getDivDiffCrafted(np.array(weeks), np.array(weights))
        print('Divided difference: {}'.format(dd[-1, -1]))
        df = rob.DataFrame({'weeks': weeks, 'weights': weights})
        rob.globalenv["y"] = df.rx2('weights')
        rob.globalenv["week"] = df.rx2('weeks')
        rob.globalenv["weight"] = rob.r('log((35-y)/y)')
        lreg = stats.lm("weight ~ week")
        a1 = lreg.rx2('coefficients')[0]
        a2 = lreg.rx2('coefficients')[1]
        # print(a1, a2)
        t = np.array(rob.r('seq(1, length(week), 0.01)'))
        z = 35 / (1 + np.exp(a1 + t * a2))
        _, ax = plt.subplots(nrows=1, ncols=1, num=1)
        plt.plot(t, z, color='purple', linewidth=3, alpha=0.65)
        plt.plot(np.array(weeks), np.array(weights), color='grey', linestyle='none', marker='o', alpha=0.7)
        plt.xlabel("Week", fontsize=12)
        plt.ylabel("Weight", fontsize=12)
        plt.annotate(
            "logistic function: " + r'$ weight = \frac{35}{1+exp({2.35-0.602*week})}$',
            fontsize=12,
            xy=(0.47, 0.48),
            xycoords='figure fraction',
        )
        plt.xticks(np.arange(10) + 1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        axe = plt.gcf()
        axe.set_size_inches(w=9.5, h=6.5)
        plt.show()
        print(t)

    def calcDivDiffs(self, file_path, file_name):
        to_dos = []
        for i in range(len(file_name)):
            f = np.transpose(pd.read_csv(file_path + file_name[i] + '.csv', sep=',', header=0))
            to_dos.append(f)
            # print(f)
            # print(pd.DataFrame(f))
        div_diffs = np.zeros([10, 3])
        for i in range(len(file_name)):
            for j in range(to_dos[i].shape[1]):
                dd = newtonDividedDifference().getDivDiffDiag(np.arange(10) + 1, to_dos[i][j])
                div_diffs[j][i] = dd[9]
                print(dd[9])
        pd.DataFrame(div_diffs).to_csv('./calc_dd.txt', sep='\t', header=None, index=False)
        return pd.DataFrame(div_diffs)

    def performTTest(self, file_path, file_name):
        to_do = pd.read_table(file_path + file_name + '.txt', sep='\t', header=None)
        for i in range(to_do.shape[1]):
            for j in range(i+1, to_do.shape[1]):
                t_stat, p_val = statisticalTest().tTestScipy(np.array(to_do[i]), np.array(to_do[j]))
                print(i+1, j+1)
                print(t_stat, p_val)

    def performFTest(self, file_path, file_name):
        to_do = pd.read_table(file_path + file_name + '.txt', sep='\t', header=None)
        t_stat, p_val = statisticalTest().fTestScipy(np.array(to_do[0]), np.array(to_do[1]), np.array(to_do[2]))
        print(t_stat, p_val)

    def poplarBCCase(self, file_path, file_name):
        to_do = np.transpose(pd.read_table(file_path + file_name + '.csv', sep=',', header=None))
        # print(to_do)
        x = np.arange(11) + 1
        # print(x)
        div_diff = []
        for i in range(1, to_do.shape[1]):
            dd = newtonDividedDifference().getDivDiffDiag(x, np.array(to_do[i].astype(np.float)))
            print(i)
            div_diff.append(dd[-1])
            print(dd[-1])
        print(div_diff)

    def poplarBCPlot(self, file_path, file_name):
        to_do = np.transpose(pd.read_table(file_path + file_name + '.csv', sep=',', header=None))
        x = np.arange(11) + 1
        fig, ax = plt.subplots()
        for i in range(1, to_do.shape[1]):
            plt.plot(x, np.array(to_do[i].astype(np.float)), marker='', color='grey', label='Populus deltoides (Poplar)' if i == to_do.shape[1] - 1 else '', linewidth=1, alpha=0.4)
        ax.set_xticks(np.arange(11) + 1)
        plt.plot(x, np.array(to_do[9].astype(np.float)), marker='', color='orange', label='Representative sample', linewidth=3.5, alpha=0.9)
        plt.legend(fontsize=10)
        plt.xlabel('Year')
        plt.ylabel('Stem Diameter(m)')
        axe = plt.gcf()
        axe.set_size_inches(w=8, h=5.5)
        plt.show()
        
    def eindexFail(self, ):
        x = np.arange(0, 10, 0.5)
        y1 = curve().logistic(x, a=1, b=1, c=1, d=5)
        y2 = curve().logistic(x, a=1, b=1, c=2.5, d=5)
        print(len(x))
        print(len(y1))
        print(len(y2))
        dd1 = newtonDividedDifference().getDivDiffDiag(x, np.array(y1))
        dd2 = newtonDividedDifference().getDivDiffDiag(x, np.array(y2))
        print(dd1[-1], dd2[-1])
        y_line = 0.1 * x
        labels = ['logistic curve 1', 'logistic curve 2', 'y=0.1*x']
        plt.plot(x, y1, marker='', color='royalblue', label=labels[0], linewidth=3, alpha=0.5)
        plt.plot(x, y2, marker='', color='r', label=labels[1], linewidth=3, alpha=0.5)
        plt.plot(x, y_line, marker='', color='grey', label=labels[2], linewidth=2, alpha=0.4)
        plt.annotate(
            "Three curves meet \nat point (0.5,0.5).",
            fontsize=9,
            xy=(0.52, 0.48),
            xycoords='figure fraction',
            xytext=(+0.52, +0.36),
            bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
            arrowprops=dict(
                arrowstyle="fancy",
                color="0.5",
                shrinkB=5,
                connectionstyle="arc3, rad=0.3",
            )
        )
        plt.legend(fontsize=10)
        plt.xlabel('Time')
        plt.ylabel('Biomass')
        axe = plt.gcf()
        axe.set_size_inches(w=8, h=5.5)
        plt.show()

    def LymphoidCurve(self, ):
        x = np.arange(0, 10.2, 0.5)
        y = np.exp(x) / (1.5 + 1 * np.exp(1.15 * (x - 5)))
        dd = newtonDividedDifference().getDivDiffDiag(x, np.array(y))
        print(dd[-1])
        fig, ax = plt.subplots(nrows=1, ncols=1, num=1)
        plt.plot(x, y, marker='', color='grey', label='Lymphoid-typed Tissue', linewidth=3, alpha=0.9)
        plt.hlines(y=68, xmin=0, xmax=10, colors="purple", linestyles="dashed", linewidth=1.5, alpha=0.6)
        plt.legend(fontsize=10)
        gca = plt.gca()
        gca.axes.get_xaxis().set_ticks([])
        gca.axes.get_yaxis().set_ticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xlabel('Time')
        plt.ylabel('Growth')
        axe = plt.gcf()
        axe.set_size_inches(w=8, h=6.5)
        plt.show()

    def ddHeatmap(self, file_path, file_name):
        to_do = np.transpose(pd.read_table(file_path + file_name + '.txt', sep='\s+', header=None))
        fig, ax = plt.subplots()
        sns.heatmap(to_do, cmap='BuPu')
        y_index = ['Group 1', 'Group 3', 'Group 3']
        ax.set_xticklabels(np.arange(10) + 1)
        ax.set_yticklabels(y_index)
        plt.xlabel('Individual')
        plt.ylabel('Group (Early QTLs Pattern)')
        axe = plt.gcf()
        axe.set_size_inches(w=8, h=5.5)
        plt.show()

    def normality(self, file_path, file_name2):
        # x_dd = [1.0333994708994844e-05, 7.10978835978841e-05, 5.643738977072275e-05, -1.5652557319224223e-05,
        #            1.6837522045855344e-05, -7.220017636684433e-06, -6.60548941798942e-05, -2.9954805996472002e-05,
        #            2.4994488536155403e-05, 3.30412257495596e-05, 7.0216049382716e-05, 3.508046737213475e-05,
        #            4.0784832451499044e-05, 5.266203703703726e-05, -8.101851851851837e-06, 6.0736331569664665e-05,
        #            4.681988536155195e-05, -1.4991181657848153e-05, -8.325066137566169e-05, 4.883156966490357e-05,
        #            4.938271604938276e-05, -8.107363315696618e-05, -9.31437389770754e-06, 2.039241622574944e-05,
        #            7.727072310405628e-05, 3.51906966490298e-05, -1.3695987654320616e-05, -1.4329805996473019e-06,
        #            3.0065035273368825e-05, 4.103284832451554e-05, -1.8022486772486615e-05, 3.3399470899470333e-05,
        #            1.1849647266309047e-06, -1.1078042328042389e-05, 1.581790123456918e-05, -7.798721340387265e-06,
        #            3.8855820105820026e-05, 1.9565696649024553e-06, 2.0392416225748628e-05, 2.9431216931217005e-05,
        #            3.7202380952381145e-05, 5.5362654320988114e-05, 9.138007054673735e-05, 5.4370590828925155e-05,
        #            -1.0141093474426953e-05, 7.864858906525519e-05, 1.4329805996472273e-05, -2.6372354497353602e-05,
        #            -4.6847442680776266e-05, -4.042658730158761e-05, 4.219025573192175e-05, -1.2400793650792143e-06,
        #            1.2676366843034504e-05, 4.309964726631406e-05, 5.798059964726615e-05, 3.0836640211640064e-05,
        #            4.819775132275137e-05, -1.5239197530863988e-05, 3.392305996472746e-05, -5.153218694885388e-05,
        #            3.4171075837743144e-05, 4.106040564373815e-05, 5.660273368606676e-05, 7.090498236331539e-05,
        #            6.429122574955892e-05, 3.5769400352733064e-05, -2.083333333333326e-05, 1.7388668430335405e-05,
        #            6.148037918871275e-05, -2.700617283949797e-06, 6.980268959435579e-05, 2.0171957671956853e-05,
        #            2.1852954144620782e-05, 3.130511463844769e-05, -8.258928571428512e-05, 5.577601410934716e-05,
        #            4.783950617283963e-05, -3.0092592592592438e-05]
        x_dd = [24.2,
34.2,
22.2,
28.5,
24.7,
30.6,
28.4,
27.7,
25.1,
22.9,
37,
19.6,
17.9,
23.2,
14.2,
34.5,
33.9,
27.7,
12.8,
33.9,
22.9,
27.3,
8.8,
41.5,
40.9,
35.4,
30.7,
16.7,
22,
22,
25.9,
38.7,
14.5,
31.2,
20.7,
12.6,
17.9,
26.1,
22.9,
34.5,
35.5,
40.2,
23.9,
22.3,
27.7,
37.4,
33.9,
32.1,
16.1,
10.7,
20.4,
14.1,
37.5,
35.7,
39,
22.6,
23.2,
14.2,
25.1,
38.3,
20.1,
23.9,
17.5,
37.7,
37,
36.1,
14.2,
35.4,
35.1,
12.3,
22,
36.7,
23.2,
22.6,
14.5,
40.9,
22.6,
34.7]
        miu1 = np.mean(x_dd)
        sigma1 = np.std(x_dd)
        print(miu1, sigma1)
        norm1 = statisticalTest().normalTest(x_dd, 'normaltest')
        print('norm1 {}'.format(norm1))
        pdf1 = probabilityDistribution().univarNorm(x_dd, miu1, sigma1)
        bind1 = pd.DataFrame({'x': x_dd, 'y': np.array(pdf1)})
        bind1 = bind1.sort_values(by=['x'], ascending=True)
        plt.plot(bind1['x'], bind1['y'], color='black', label='Probability of the NT-based phenotypic data', linestyle='solid', alpha=0.6)
        plt.hist(x_dd, bins=30, density='norm', label='NT-based phenotypic data', facecolor='skyblue', ec="grey", alpha=0.6)
        plt.text(-2.3, 0.3, r'$\mu=3.42e-17,\ \sigma=1$')
        plt.legend(fontsize=10)
        axe = plt.gcf()
        axe.set_size_inches(w=8, h=6.5)
        plt.show()

        x_ori = pd.read_table(file_path + file_name2 + '.csv', sep=',', header=0)
        x_ori_ = np.array(x_ori['11th Year'])
        miu2 = np.mean(x_ori_)
        sigma2 = np.std(x_ori_)
        print(miu2, sigma2)
        norm2 = statisticalTest().normalTest(x_ori_, 'normaltest')
        print('norm1 {}'.format(norm2))
        pdf2 = probabilityDistribution().univarNorm(x_ori_, miu2, sigma2)
        bind2 = pd.DataFrame({'x': x_ori_, 'y': np.array(pdf2)})
        bind2 = bind2.sort_values(by=['x'], ascending=True)
        plt.plot(bind2['x'], bind2['y'], color='black', label='Probability of the original phenotypic data', linestyle='solid', alpha=0.6)
        plt.hist(x_ori_, bins=30, density='norm', label='Original phenotypic data', facecolor='skyblue', ec="grey", alpha=0.6)
        plt.text(1.7, 0.28, r'$\mu=-4.73e-16,\ \sigma=1$')
        plt.legend(fontsize=10)
        axe = plt.gcf()
        axe.set_size_inches(w=8, h=6.5)
        plt.show()

    def fourInterpos(self, ):
        x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        y = [0, 1, 1.5, 2, 6, 10.5, 11, 15, 31, 33, 34]
        x_new = np.linspace(0, 10, 100)
        # fig, ax = plt.subplot(nrows=1, ncols=3, num=1)
        plt.scatter(x, y, marker='*')
        y1 = interpolation().linear(x, y)
        y2, _ = interpolation().csi(x, y)
        y3, _ = interpolation().pchip(x, y)
        y_new1 = y1(x_new)
        y_new2 = y2(x_new)
        y_new3 = y3(x_new)
        plt.subplot(221)
        plt.scatter(x, y, marker='^', lw=3, c='cadetblue', alpha=0.6)
        plt.plot(x_new, y_new1, c='cadetblue', label='Linear', lw=3, alpha=0.6)
        plt.title('(a)')
        plt.ylabel('Growth')
        plt.legend(fontsize=12)
        plt.subplot(222)
        plt.scatter(x, y, marker='^', lw=3, c='darkorange', alpha=0.6)
        plt.plot(x_new, y_new2, c='darkorange', label='CSI', lw=3, alpha=0.6)
        plt.title('(b)')
        plt.legend(fontsize=12)
        plt.subplot(223)
        plt.scatter(x, y, marker='^', lw=3, c='darkgreen', alpha=0.6)
        plt.plot(x_new, y_new3, c='darkgreen', label='PCHIP', lw=3, alpha=0.6)
        plt.title('(c)')
        plt.xlabel('Time')
        plt.ylabel('Growth')
        plt.legend(fontsize=12)

        ### /* newton interpolation */ ###
        plt.subplot(224)
        y_new4 = interpolation().newton(x, y, x_new)
        plt.scatter(x, y, marker='^', lw=3, c='dimgray', alpha=0.6)
        plt.plot(x_new, y_new4, c='dimgray', label='newton', lw=3, alpha=0.6)
        plt.title('(d)')
        plt.xlabel('Time')
        plt.annotate(
            "Leading to runge phenomenon \nof higher-order interpolation\n at last three points.",
            fontsize=9,
            xy=(0.85, 0.40),
            xycoords='figure fraction',
            xytext=(+0.60, +0.32),
            bbox=dict(boxstyle="round4,pad=.5", fc="0.8"),
            arrowprops=dict(
                arrowstyle="fancy",
                color="0.5",
                shrinkB=5,
                connectionstyle="arc3, rad=0.3",
            ),
        )
        axe = plt.gcf()
        axe.set_size_inches(w=10.5, h=8.0)
        plt.legend(fontsize=12)
        plt.show()

    def lowMidUp(self, ):
        upper = np.zeros([10, 501])
        middle = np.zeros([10, 501])
        lower = np.zeros([10, 501])
        for i in range(10):
            for j in range(501):
                step = 0.02 * j
                upper[i][j] = curve().logistic(step, a=1, b=1, c=2+0.3*i, d=2)
                middle[i][j] = curve().logistic(step, a=1, b=1, c=2.3+0.3*i, d=5)
                lower[i][j] = curve().logistic(step, a=1, b=1, c=2.3+0.3*i, d=8)
        x = np.arange(0, 10.02, 0.02)
        _, ax = plt.subplots(nrows=1, ncols=1, num=1)
        for i in range(10):
            plt.plot(x, upper[i], c='darkgreen', label='Group 1' if i == 9 else '', lw=1.5, alpha=0.6)
            plt.plot(x, middle[i], c='cadetblue', label='Group 2' if i == 9 else '', lw=1.5, alpha=0.6)
            plt.plot(x, lower[i], c='dimgray', label='Group 3' if i == 9 else '', lw=1.5, alpha=0.9)
        plt.legend(fontsize=10)
        plt.xlabel('Time', fontsize=12)
        plt.ylabel('Growth Degree', fontsize=12)
        T = np.arange(0, 1.2, 0.2)
        print(T)
        ax.set_yticks(T)
        y_index = ['0', '20%', '40%', '60%', '80%', '100%']
        ax.set_yticklabels(y_index)
        axe = plt.gcf()
        axe.set_size_inches(w=11, h=6.5)
        plt.show()

        upper_sample = np.zeros([10, 10])
        mid_sample = np.zeros([10, 10])
        lower_sample = np.zeros([10, 10])
        for i in range(10):
            c = 0
            for j in range(10):
                c = c + 50
                upper_sample[i][j] = upper[i][c]
                mid_sample[i][j] = middle[i][c]
                lower_sample[i][j] = lower[i][c]
        return upper_sample


if __name__ == "__main__":
    p = ntBoard()

    # print(p.mouseCaseStudy())

    # print(p.calcDivDiffs(
    #     file_path='./data/simu/',
    #     file_name=['upper_points', 'mean_points', 'lower_points']
    # ))

    # print(p.performTTest(
    #     file_path='./',
    #     file_name='calc_dd'
    # ))

    # print(p.performFTest(
    #     file_path='./',
    #     file_name='calc_dd'
    # ))

    # print(p.poplarBCCase(
    #     file_path='./',
    #     file_name='populus.BC.pheno'
    # ))

    # print(p.poplarBCPlot(
    #     file_path='./',
    #     file_name='populus.BC.pheno'
    # ))

    # print(p.eindexFail())

    # print(p.LymphoidCurve())

    # print(p.ddHeatmap(
    #     file_path='./',
    #     file_name='calc_dd'
    # ))

    print(p.normality(
        file_path='./',
        file_name2='populus.BC.pheno'
    ))

    # print(p.fourInterpos())

    # print(p.lowMidUp())
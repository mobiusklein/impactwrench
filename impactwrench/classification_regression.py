import numpy as np
import seaborn as sns
from mlxtend.plotting import plot_decision_regions
from sklearn.svm import SVR
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.preprocessing import PolynomialFeatures
from sklearn.svm import SVC
from sklearn.svm import LinearSVR
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn import tree
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from mlxtend.plotting import plot_decision_regions
import matplotlib.gridspec as gridspec
import itertools
from sklearn.externals.six import StringIO  
from IPython.display import Image  
from sklearn.tree import export_graphviz
import pydotplus

def main():

       def linear_svm(X_train, y_train):
            svm_clf = Pipeline((( "scaler", StandardScaler()), 
            ("linear_svc", LinearSVC(C = 1, loss = "hinge"))))
            return svm_clf.fit(X, y)

#for non linearly seperable datasets, use the following

      def non_linear_svm(Xtrain, y_train):       
            polynomial_svm_clf = Pipeline((("poly_features", PolynomialFeatures(degree = 3)),
            ("scaler", StandardScaler()),
            ("svm_clf", LinearSVC(C = 10, loss = "hinge")))
            return polynomial_svm_clf.fit(X_train, y_train)


#kernel trick for using linearSVC on non-linear data

        def polykernel_svm(X_train, y_train):
            poly_kernel_svm_clf = Pipeline((("scaler", StandardScaler()),
            ("svm_clf", SVC(kernel = "poly", degree = 3, coef0=1, C = 5))))
            return poly_kernel_svm_clf.fit(X,y)

    #using radial bias functions
        def rbfkernel(X_train, y_train):
            rbf_kernel_svm_clf = Pipeline((("scaler", StandardScaler()),
                ("svm_clf", SVC(kernel = "rbf", gamma = 5, C = 0.001)) ))
            return rbf_kernel_svm_clf.fit(X,y)

#linear regression

# multiple output regressor

        def multiple_output_RF_regressor(X_train, y_train):
        "This function does a simple multioutput Random Forest Regression, to visualize the output of the function"
        "refer to the plots.py which contains all the plotting procedures."
        
        max_depth = [10, 30]
        random_state= [0, 1, 2]

        for i in max_depth:
                 regr_multirf = MultiOutputRegressor(RandomForestRegressor(max_depth=i,random_state=0))
          return  regr_multirf.fit(X_train, y_train)
        
        for i in random_state:
             regr_rf = RandomForestRegressor(max_depth=max_depth, random_state=i)
            return regr_rf.fit(X_train, y_train)

        # Predict on new data
        y_multirf = regr_multirf.predict(X_test)
        y_rf = regr_rf.predict(X_test)

        plt.figure()
        s = 50
        a = 0.4
        plt.scatter(y_test[:, 0], y_test[:, 0], edgecolor='k',
                    c="navy", s=s, marker="s", alpha=a, label="Data")
        plt.scatter(y_multirf[:, 0], y_multirf[:, 0], edgecolor='k',
                    c="cornflowerblue", s=s, alpha=a,
                    label="Multi RF score=%.2f" % regr_multirf.score(X_test, y_test))
        plt.scatter(y_rf[:, 0], y_rf[:, 0], edgecolor='k',
                    c="c", s=s, marker="^", alpha=a,
                    label="RF score=%.2f" % regr_rf.score(X_test, y_test))
        plt.xlim([-6, 6])
        plt.ylim([-6, 6])
        plt.xlabel("target 1")
        plt.ylabel("target 2")
        plt.title("Comparing random forests and the multi-output meta estimator")
        plt.legend()
        plt.show()

        #simple linear regression
        svm_reg= LinearSVR(epsilon = 1.5)
        svm_reg.fit(X,y)


#plotting decision boundary between classes
       def one_dimensional_decision_boundary(X_train, y_train):
           svr = SVR(C=0.5, kernel='linear')
           svr.fit(X, y)
           plot_decision_regions(X, y, clf=svr, legend=2)

            plot_decision_regions(X, y, clf=svr, legend=2)

            # Adding axes annotations
            plt.xlabel(' Intensity [units]')
            plt.title('SVR')

            plt.show()

# multiple classifiers
       def classifier_comparison(X_train, y_train):
            clf1 = LogisticRegression(random_state=1)
            clf2 = RandomForestClassifier(random_state=1)
            clf3 = GaussianNB()
            clf4 = SVC()    

            gs = gridspec.GridSpec(2, 2)

            fig = plt.figure(figsize=(10,8))

            labels = ['Logistic Regression', 'Random Forest', 'Naive Bayes', 'SVM']
            for clf, lab, grd in zip([clf1, clf2, clf3, clf4],
                                    labels,
                                    itertools.product([0, 1], repeat=2)):

                clf.fit(X, y)
                ax = plt.subplot(gs[grd[0], grd[1]])
                fig = plot_decision_regions(X=X, y=y, clf=clf, legend=2)
                plt.title(lab)

            plt.show()

# grid of decision slices

       def decision_slices(X_train, y_train):
            svm = SVC()
            svm.fit(X, y)

            # Plotting decision regions
            fig, axarr = plt.subplots(2, 2, figsize=(10,8), sharex=True, sharey=True)
            values = [ 500, 1000, 2000]
            width = 500
            for value, ax in zip(values, axarr.flat):
                plot_decision_regions(X, y, clf=svm,
                                    filler_feature_values={2: value},
                                    filler_feature_ranges={2: width},
                                    legend=2, ax=ax)
                ax.set_xlabel('Feature 1')
                ax.set_ylabel('Feature 2')
                ax.set_title('Feature 3 = {}'.format(value))

            # Adding axes annotations
            fig.suptitle('SVM on Etanercept')
            plt.show()

  # multiple learners 

        dict_classifiers = {
            "Logistic Regression": LogisticRegression(),
            "Nearest Neighbors": KNeighborsClassifier(),
            "Linear SVM": SVC(),
            "Gradient Boosting Classifier": GradientBoostingClassifier(n_estimators=1000),
            "Decision Tree": tree.DecisionTreeClassifier(),
            "Random Forest": RandomForestClassifier(n_estimators=1000),
            "Neural Net": MLPClassifier(alpha = 1),
            "Naive Bayes": GaussianNB(),
            "AdaBoost": AdaBoostClassifier(),
            "QDA": QuadraticDiscriminantAnalysis(),
            "Gaussian Process": GaussianProcessClassifier() 
            }  

        import time
        def batch_classify(X_train, Y_train, X_test, Y_test, no_classifiers = 5, verbose = True):
        
        dict_models = {}
        for classifier_name, classifier in list(dict_classifiers.items())[:no_classifiers]:
            t_start = time.clock()
            classifier.fit(X_train, Y_train)
            t_end = time.clock()
        
            t_diff = t_end - t_start
            train_score = classifier.score(X_train, Y_train)
            test_score = classifier.score(X_test, Y_test)
            
            dict_models[classifier_name] = {'model': classifier, 'train_score': train_score, 'test_score': test_score, 'train_time': t_diff}
            if verbose:
                print("trained {c} in {f:.2f} s".format(c=classifier_name, f=t_diff))
        return dict_models
 
 
 
        def display_dict_models(dict_models, sort_by='test_score'):
            cls = [key for key in dict_models.keys()]
            test_s = [dict_models[key]['test_score'] for key in cls]
            training_s = [dict_models[key]['train_score'] for key in cls]
            training_t = [dict_models[key]['train_time'] for key in cls]
            
            df_ = pd.DataFrame(data=np.zeros(shape=(len(cls),4)), columns = ['classifier', 'train_score', 'test_score', 'train_time'])
            for ii in range(0,len(cls)):
                df_.loc[ii, 'classifier'] = cls[ii]
                df_.loc[ii, 'train_score'] = training_s[ii]
                df_.loc[ii, 'test_score'] = test_s[ii]
                df_.loc[ii, 'train_time'] = training_t[ii]
            
            display(df_.sort_values(by=sort_by, ascending=False))

# classification trees

         def classification_tree(X_train, y_train):
             tree_clf = DecisionTreeClassifier(max_depth = 2)
             tree_clf.fit(X_train,Y_train)

             from sklearn.tree import export_graphviz
             export_graphviz(tree_clf,out_file = "classification_tree.dot",max_depth = None,feature_names = None,
                rounded = True,
                filled = True)


# ensemble learners

        def ensemble_learner(X_train, y_train):
            log_clf = LogisticRegression()
            rnd_reg = RandomForestRegressor()
            svm_reg = LinearSVR()

            voting_clf = VotingClassifier(
            estimators = [('lr', log_clf), ('rf', rnd_reg), ('svc', svm_reg)],
            voting = 'soft')  

#bootstrap aggegation
        from sklearn.ensemble import BaggingRegressor
        from sklearn.tree import DecisionTreeRegressor
        from sklearn.metrics import accuracy_score

       def bagging_regressor(X_train, y_train):
           for max_number_samples = len(X_train):
                if n_estimators  < len(y_train):
               bag_reg = BaggingRegressor(DecisionTreeRegressor(), n_estimators,
               max_samples = 50, bootstrap = True, n_jobs = -1)
               bag_reg.fit(X_train, Y_train)    # n_jobs denotes the number of processors dedicated to the the process
               y_pred = bag_reg.predict(X_test)
            else: # out of bag evaluation
                bag_reg = BaggingRegressor(
                DecisionTreeRegressor(), n_estimators,
                bootstrap = True, n_jobs = -1, oob_score = True)
                bag_reg.fit(X_train, Y_train)
                bag_reg.oob_score_

#seqeuencial decision tree regressor/ gradient boosting

        def decision_tree_stacking(X_train, y_train):
            " a number of first level individual learners are generated from training data set" 
            "and then combined by a metalearner"
                tree_reg1 = DecisionTreeRegressor(max_depth = 2)
                tree_reg1.fit(X_train, Y_train)

                # now train a second DecisionTreeRegressor on the residual errors made by the first predictor

                y2 = Y_train - tree_reg1.predict(X_train)
                tree_reg2 = DecisionTreeRegressor(max_depth = 2)
                tree_reg2.fit(X_train, y2)

                # now train a third regressor  on the residual errors made by the second predictor

                y3 = y2-tree_reg.predict(X_train)
                tree_reg3 = DecisionTreeRegressor(max_depth = 2)
                tree_reg3.fit(X_train,y3)

                #the ensemble of three trees can now be used to make new predictions

                y_pred = sum(tree.predict(X_test) for tree in (tree_reg1, tree_reg2, tree_reg3))

# using the early stopping algorithm as a regularization method

       def early_stopping(X_train, y_train):
            gbrt = GradientBoostingRegressor(max_depth = 2, n_estimators = 120)
            gbrt.fit(X_train, y_train)

            errors = [mean_squared_error(y_val, y_pred)
                    for y_pred in gbrt.staged_predict(X_val)]
            bst_n_estimators = np.argmin(errors)

            gbrt_best = GradientBoostingRegressor(max_depth = 2, n_estimators= bst_n_estimators)
            return gbrt_best.fit(X_train, y_train)

# getting validation error

        def validation_error(X_train, y_train):
            gbrt = GradientBoostingRegressor(max_depth = 2, warm_start = True)

            min_val_error = float("inf")
            error_going_up = 0
            for n_estimators in range(1, 120):
                gbrt.n_estimators = n_estimators
                gbrt.fit(X_train, y_train)
                y_pred = gbrt.predict(X_val)
                val_error = mean_squared_error(y_val, y_pred)
                if val_error < min_val_error:
                    min_val_error = val_error
                    error_going_up = 0
                else:
                    error_going_up += 1
                    if error_going_up ==5:
                        break #early stopping

#plotting a decision surface

  # Parameters
        def plot_decision_surface(X_train, y_train):
            n_classes = 3
            plot_colors = "ryb"
            plot_step = 0.02

            for pairidx, pair in enumerate([[0, 1], [0, 2], [0, 3],
                                            [1, 2], [1, 3], [2, 3]]):
                

                # Train
                clf = DecisionTreeClassifier().fit(X_train, y_train)

                # Plot the decision boundary
                plt.subplot(2, 3, pairidx + 1)

                x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
                y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
                xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step),
                                    np.arange(y_min, y_max, plot_step))
                plt.tight_layout(h_pad=0.5, w_pad=0.5, pad=2.5)

                Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
                Z = Z.reshape(xx.shape)
                cs = plt.contourf(xx, yy, Z, cmap=plt.cm.RdYlBu)

                plt.xlabel('expected m/Z')
                plt.ylabel('pvalues')

                # Plot the training points
                for i, color in zip(range(n_classes), plot_colors):
                    idx = np.where(y == i)
                    plt.scatter(X[idx, 0], X[idx, 1], c=color,
                                cmap=plt.cm.RdYlBu, edgecolor='black', s=15)

            plt.suptitle("Decision surface of a decision tree using paired features")
            #plt.legend(loc='lower right', borderpad=0, handletextpad=0)
            plt.axis("tight")
            plt.show()
    
if __name__ == '__main__':
    main()   
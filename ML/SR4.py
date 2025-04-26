import numpy as np
from sklearn.model_selection import KFold
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import learning_curve, ParameterGrid, train_test_split, GridSearchCV


# To automate labeling graphs, and load data from text files
# Change file name depending on property or descriptor being tested
directory = "Combined Testing finalization/"
descriptor_filename = ["DescriptorV1.txt", "DescriptorV1.txt", "DescriptorEditedElectroAff.txt", "DescriptorEditedElectroNeg.txt"]
z_desc_filename = ["AtomicNumber.txt", "AtomicNumber.txt", "ZEditedElectroAff.txt","ZEditedElectroNeg.txt"]
property_filename = ["AtomicRadius.txt", "FirstIonEn.txt", "EditedElectroAff.txt", "EditedElectroNeg.txt"]


# Making parameter ranges
alphas = [10**i for i in range(-10, 8)]
gammas = [10**i for i in range(-10, 8)]

# Checking parameters for each type of ridge regression.
param_grid = [
    {'alpha': alphas,
     'kernel': ['linear']},
    {'alpha': alphas, 'gamma': gammas,
     'kernel': ['rbf', 'poly']}
]

def descriptor_testing(X, y):
    # Splits into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=.80, shuffle=True, random_state=1) 
    # Cross validation
    search = GridSearchCV(KernelRidge(), param_grid, n_jobs=-1, cv=5).fit(X_train,y_train)
    # Prediction
    best_model = search.best_estimator_
    y_pred_all = best_model.predict(X)

    # Prediction Performance of the entire model
    mae_test = mean_absolute_error(y, y_pred_all)
    r2_test = r2_score(y, y_pred_all)
    mse_test = mean_squared_error(y, y_pred_all)

    return y_pred_all, mae_test, r2_test, mse_test

def Z_testing(X, y):
    #Splits into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=.80, shuffle=True, random_state=1) 
    X_train = X_train.reshape(-1,1)
    X_test = X_test.reshape(-1,1)
    # Cross validation
    search = GridSearchCV(KernelRidge(), param_grid, n_jobs=-1, cv=5).fit(X_train,y_train)
    # Prediction
    best_model = search.best_estimator_
    y_pred_all = best_model.predict(X.reshape(-1,1))

    # Prediction Performance of the entire model
    mae_test = mean_absolute_error(y, y_pred_all)
    r2_test = r2_score(y, y_pred_all)
    mse_test = mean_squared_error(y, y_pred_all)

    return y_pred_all, mae_test, r2_test, mse_test


with open(f"{directory} 80% test Descriptor and Z testing.txt", "w") as f:
    for index, names in enumerate(descriptor_filename):
        X = np.loadtxt(descriptor_filename[index])
        Z = np.loadtxt(z_desc_filename[index])
        y = np.loadtxt(property_filename[index])
        X_filename = descriptor_filename[index].removesuffix('.txt')
        y_filename = property_filename[index].removesuffix('.txt')

        # Assigning values from descriptor and z testing
        d_test = descriptor_testing(X, y)
        d_pred = d_test[0]
        d_mae = d_test[1]
        d_r2 = d_test[2]
        d_mse = d_test[3]

        z_test = Z_testing(Z,y)
        z_pred = z_test[0]
        z_mae = z_test[1]
        z_r2 = z_test[2]
        z_mse = z_test[3]

        # Cases for presenting results
        match property_filename[index]:
            case "EditedElectroNeg.txt":
                prop_test = "Electronegativity"
            case "FirstIonEn.txt":
                prop_test = "First Ionization Energy"
            case "EditedElectroAff.txt":
                prop_test = "Electron Affinity"
            case "AtomicRadius.txt":
                prop_test = "vdW Atomic Radius"
        
        match property_filename[index]:
            case "EditedElectroNeg.txt":
                units = ""
            case "FirstIonEn.txt":
                units = "(eV)"
            case "EditedElectroAff.txt":
                units = "(eV)"
            case "AtomicRadius.txt":
                units = "(pm)"

        # Saving data
        f.write(f"Results from 80 {prop_test} testing:\n")
        f.write("------------------\n")
        f.write(f"Mean absolute error: {d_mae}\n")
        f.write(f"Mean square error: {d_mse}\n")
        f.write(f"R^2 value: {d_r2}\n")
        f.write("------------------\n")
        f.write("Results from Z testing:\n")
        f.write("------------------\n")
        f.write(f"Mean absolute error: {z_mae}\n")
        f.write(f"Mean square error: {z_mse}\n")
        f.write(f"R^2 value: {z_r2}\n")
        f.write("------------------\n")

        # Plotting predicted vs actual
        fig, axes = plt.subplots(1, 2, figsize=(12, 6)) # figsize=? can be used to adjust size
        ideal_line = np.linspace(min(y), max(y))
        d_resid = y - d_pred
        z_resid = y - z_pred
        std_residual_d = np.std(d_resid)
        upper_bound = ideal_line + std_residual_d
        lower_bound = ideal_line - std_residual_d

        # Using Z
        axes[0].scatter(y, z_pred, color='blue', label='Z', s=9)
        axes[0].set_xlabel(f"Actual {prop_test} {units}", fontsize=14)
        axes[0].set_ylabel(f"Estimated {prop_test} {units}", fontsize=14)
        axes[0].plot([min(y), max(y)], [min(y), max(y)], 'r--', label='Ideal')
        axes[0].plot(ideal_line, upper_bound, color='gray', linestyle='--', linewidth=0.7, label='+/- 1 std from Ideal')
        axes[0].plot(ideal_line, lower_bound, color='gray', linestyle='--', linewidth=0.7)
        axes[0].legend(loc="best")
        axes[0].text(0.05, 0.95, f'MAE: {z_mae:.2f}', transform=axes[0].transAxes, verticalalignment="top", fontsize=12)
        axes[0].set_title("Using Z", fontsize=14)
        # Using Descriptor
        axes[1].scatter(y, d_pred, color='black', label='RCDD', s=9)
        axes[1].set_xlabel(f"Actual {prop_test} {units}", fontsize=14)
        axes[1].set_ylabel(f"Estimated {prop_test} {units}", fontsize=14)
        axes[1].plot([min(y), max(y)], [min(y), max(y)], 'r--', label='Ideal')
        axes[1].plot(ideal_line, upper_bound, color='gray', linestyle='--', linewidth=0.7, label='+/- 1 std from Ideal')
        axes[1].plot(ideal_line, lower_bound, color='gray', linestyle='--', linewidth=0.7)
        axes[1].legend(loc="best")
        axes[1].text(0.05, 0.95, f'MAE: {d_mae:.2f}', transform=axes[1].transAxes, verticalalignment="top", fontsize=12)
        axes[1].set_title("Using RCDD", fontsize=14)

        # General plot items
        plt.suptitle(f"Correlation of predictions to {prop_test} data", fontsize=14, fontweight='bold')
        plt.tight_layout
        plt.savefig(f"{directory}SR4 80 Correlation of predictions to {prop_test} data")
        plt.show()

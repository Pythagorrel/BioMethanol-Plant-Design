module Final_Margules_Fit

export perform_margules_fit

# 1. Load Modules

using Reexport

@reexport using VLE_MeH2O

@reexport using GeRT   

@reexport using GeRTx1x2



using Printf


# --- Your Curve Fitting Function of Choice---

function reglin(x,y) 

    n=length(x)

    length(y)==n||error("x and y must be the same length")

    sx=sum(x);sy=sum(y)

    sx2=sum(x.^2);sy2=sum(y.^2);sxy=sum(x.*y)

    mSlo=(n .*sxy .- sx .*sy)/(n .*sx2 .-sx.^2) #gives gradient

    clnt=sy ./n .-mSlo .*sx ./n #gives intercept

    return mSlo, clnt

end


# --- Main Execution Logic ---

function perform_margules_fit()

    

    # 2. Get All Data

    x_all = VLE_MeH2O.x_methanol

    y_all = GeRTx1x2.get_GeRTx1x2()

    

    # 3. Slice Data to Remove Endpoints

    # We ignore the first and last points (x=0 and x=1) 

    # to avoid the "division by epsilon" spikes.

    x_fit = x_all[2:end-1]

    y_fit = y_all[2:end-1]

    

    # 4. Perform Linear Regression

    slope, intercept = reglin(x_fit, y_fit)

    

    # 5. Calculate Margules Parameters

    # Theory: Ge/(RTx1x2) = (A21 - A12)x1 + A12

    # Therefore:

    # Intercept = A12

    # Slope = A21 - A12  =>  A21 = Slope + Intercept

    

    A12 = intercept

    A21 = slope + intercept

    

    #= 6. Output Results

    println("Margules Model Fitting Results")

    println("------------------------------")

    @printf("Slope (m):      %6.4f\n", slope)

    @printf("Intercept (c):  %6.4f\n", intercept)

    println("------------------------------")

    println("Calculated Binary Parameters:")

    @printf("A12 = %6.4f  (Limit as x1 -> 0)\n", A12)

    @printf("A21 = %6.4f  (Limit as x1 -> 1)\n", A21)=#

    

    return A12, A21

end

end
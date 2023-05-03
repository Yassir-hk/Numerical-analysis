/**************************************************************************************************
 *                                     Handlers and utilities                                     *
 **************************************************************************************************/

const regexVar = /^[a-zA-Z]+$/;
const regexNum = /^[\d.]+$/;
const regexOpt = /^[+-]$/;

/**
 * Function to change the number of equations.
 * @param {Number | Boolean} increment
 */

function changeEquationNumber(increment) {
  const numberOfEquationsTracker = document.getElementById("equationNumber");
  const numberOfEquations = parseInt(numberOfEquationsTracker.innerHTML) + increment;

  if (numberOfEquations < 0 || numberOfEquations > 6) {
    return false;
  }

  numberOfEquationsTracker.innerHTML = numberOfEquations;
  const equationsContainer = document.getElementById("equationsContainer");

  if (increment == 1) {
    equationsContainer.appendChild(createEquation(numberOfEquations));
  }
  else {
    equationsContainer.removeChild(equationsContainer.lastChild);
  }
}

/**
 * Function to create and equation.
 * @param {Number} equationNumber 
 */

function createEquation(equationNumber) {
  const equation = document.createElement("div");
  const equationInput = document.createElement("input");
  const equationName = document.createTextNode("Equation " + equationNumber);

  equationInput.type = "text";
  equation.className = "equation";
  equationInput.className = "equationInput equation" + equationNumber;
  
  equation.appendChild(equationName);
  equation.appendChild(equationName);
  equation.appendChild(equationInput);

  return equation;
}

/**
 * Function to check if a given equation is following the correct format.
 * @param {Array | null} equation
 * @returns {Boolean}
 */

function isValidEquation(equation) {
  let lastChar = false;
  let countEquales = 0;

  for (let char of equation) {
    const isValidNum = regexNum.test(char) && (regexOpt.test(lastChar) || regexNum.test(lastChar) || lastChar == false);
    const isValidOpt = regexOpt.test(char) && regexVar.test(lastChar);
    const isValidVar = regexVar.test(char) && regexNum.test(lastChar);

    if (!isValidNum && !isValidOpt && !isValidVar) {
      return false;
    }

    countEquales += (char === '=');
    lastChar = char;
  }

  return countEquales == 1;
}

/**
 * Function to verity the equations input.
 * @returns {Array | boolean}
 */

function fetchSystemData() {
  const equations = Array.from(document.getElementsByClassName("equationInput"));
  const variablesOfSystemSet = new Set();
  const systemData = [];

  for (let equation of equations) {
    const equationString = equation.value.replace(/\s+/g, '');

    if (!isValidEquation(equationString)) {
      return false;
    }

    const systemEquation = {};
    let coefficient = '', operator = '';

    for (let char of equationString) {
      if (regexNum.test(char)) {
        coefficient += char;
      }
      else if (regexVar.test(char)) {
        if (!variablesOfSystemSet.has(char)) {
          variablesOfSystemSet.add(char);
        }
        if (!systemEquation[char]) {
          systemEquation[char] = 0;
        }
        if (coefficient.length == 0) {
          coefficient = '1';
        }

        systemEquation[char] += parseFloat(operator + coefficient);
        coefficient = '';
        operator = '';
      }
      else if (regexOpt.test(char)) {
        operator = char;
      }
    }
    if (coefficient.length) {
      systemEquation['='] = parseFloat((operator == '=' ? '' : operator) + coefficient);
    }
    systemData.push(systemEquation);
  }

  // Check if the number of variables is greater than the number of equations
  if (variablesOfSystemSet.size > equations.length) {
    return false;
  }

  // Add coefficients with value 0
  const variablesOfSystemArray = [...variablesOfSystemSet].sort();
  const matrixOfSystem = [];

  for (let systemEquation of systemData) {
    const equationCoefficients = [];

    for (let variable of variablesOfSystemArray) {
      if (variable in systemEquation) {
        equationCoefficients.push(systemEquation[variable]);
      }
      else {
        equationCoefficients.push(0);
      }
    }
    equationCoefficients.push(systemEquation['=']);
    matrixOfSystem.push(equationCoefficients);
  }

  return [
    matrixOfSystem,
    variablesOfSystemArray
  ];
}

/**
 * Function that process the input equations and use the selected method to provide the resultVector.
 * @returns {Boolean | void}
 */

function compute() {
  const systemData = fetchSystemData();
  const computationMethod = parseInt(document.getElementById("computationMethod").value);

  // Case of invalid inputs or invalid system
  if (systemData == false) {
    alert("This system cannot be processed");
    return false;
  }

  const resultVectorContainer = Array.from(document.getElementsByClassName("container-2"));
  resultVectorContainer[0].style.visibility = "visible";

  // Use the selected method
  switch (computationMethod) {
    case 1: {
      const resultVectorData = gaussElemination(systemData[0]);
      writeresultVector(resultVectorData, systemData[1]);
      break;
    }
    case 2: {
      const resultVectorData = luDecomposition(systemData[0]);
      writeresultVector(resultVectorData, systemData[1]);
      break;
    }
    case 3: {
      break;
    }
  }
}

/**
 * Function to write the solution of the system in the resultVector container.
 * @param {Array | Boolean} systemSolution
 * @param {Array} variablesOfSystem
 * @returns {Boolean}
 */

function writeresultVector(systemSolution, variablesOfSystem) {
  const resultVectorContainer = document.getElementById("resultVectorContainer");

  if (systemSolution == false) {
    resultVectorContainer.innerHTML = "This system has infinity of solutions or no solution.";
    return false;
  }

  let resultVectorContent = '';
  let currtVariable = 0;

  for (let value of systemSolution) {
    resultVectorContent += variablesOfSystem[currtVariable++] + " = " + value + "<br><br>";
  }

  resultVectorContainer.innerHTML = resultVectorContent;
}

/**************************************************************************************************
 *                                      Gauss elemination                                         *
 **************************************************************************************************/

/**
 * Function is used to solve a system of linear equations using the Gaussian elimination method. 
 * 1) Triangularizing the system
 * 2) solve the system
 * 
 * @param {Array} system 
 * @returns {Array | Boolean}
 */

function gaussElemination(system) {
  const numberOfRows = system.length;
  const numberOfCols = system[0].length;

  // Triangularization of the system
  for (let i = 0; i < numberOfRows; i++) {
    let pivotIndex = i;

    while (pivotIndex < numberOfRows && system[pivotIndex][i] == 0) {
      pivotIndex++;
    }

    if (pivotIndex < numberOfRows && pivotIndex != i) {
      [system[pivotIndex], system[i]] = [system[i], system[pivotIndex]];
    }
    else if (pivotIndex == numberOfRows) {
      continue;
    }

    for (let j = i + 1; j < numberOfRows; j++) {
      const factor = system[j][i] / system[i][i];

      for (let k = i; k < numberOfCols; k++) {
        system[j][k] -= system[i][k] * factor;
      }
    }
  }

  return upperTriangularSystemSolver(system);
}

/**
 * Function that solve the upper triangular system
 * @param {Array} system
 */

function upperTriangularSystemSolver(system) {
  const numberOfRows = system.length;
  const numberOfCols = system[0].length;

  if (numberOfRows + 1 != numberOfCols) {
    return false;
  }

  let resultVector = new Array(numberOfRows);
  
  for (let i = numberOfRows - 1; i >= 0; i--) {
    resultVector[i] = system[i][numberOfCols - 1];

    for (let j = numberOfRows - 1; j > i; j--) {
      resultVector[i] -= resultVector[j] * system[i][j];
    }

    resultVector[i] /= system[i][i];
  }

  return resultVector;
}

/**
 * Function that solve the lower triangular system
 * @param {Array} system 
 */

function lowerTriangularSystemSolver(system) {
  const numberOfRows = system.length;
  const numberOfCols = system[0].length;

  if (numberOfRows + 1 != numberOfCols) {
    return false;
  }

  let resultVector = new Array(numberOfRows);

  for (let i = 0; i < numberOfRows; i++) {
    resultVector[i] = system[i][numberOfCols - 1];

    for (let j = 0; j < i; j++) {
      resultVector[i] -= resultVector[j] * system[i][j];
    }

    resultVector[i] /= system[i][i];
  }

  return resultVector;
}

/**************************************************************************************************
 *                                      LU decomposition                                          *
 **************************************************************************************************/

/**
 * Function to get the lower and the upper matrices
 * @param {Array} system
 * @returns {Array | Boolean}
 */

function luDecomposition(system) {
  const dimension = system.length;
  const lowerTriangularMatrix = [];
  const upperTriangularMatrix = Array.from(system);

  // Initializing the lower matrix  
  for (let i = 0; i < dimension; i++) {
    const matrixRow = [];

    for (let j = 0; j < dimension; j++) {
      matrixRow.push((i == j) ? 1 : 0);
    }
    
    lowerTriangularMatrix.push(matrixRow);
  }

  for (let i = 0; i < dimension; i++) {
    let pivotIndex = i;

    while (pivotIndex < dimension && upperTriangularMatrix[pivotIndex][i] == 0) {
      pivotIndex++;
    }

    if (pivotIndex < dimension && pivotIndex != i) {
      [upperTriangularMatrix[pivotIndex], upperTriangularMatrix[i]] = [upperTriangularMatrix[i], upperTriangularMatrix[pivotIndex]];
    }
    else if (pivotIndex == dimension) {
      continue;
    }

    for (let j = i + 1; j < dimension; j++) {
      const factor = upperTriangularMatrix[j][i] / upperTriangularMatrix[i][i];

      for (let k = i; k < dimension; k++) {
        upperTriangularMatrix[j][k] -= upperTriangularMatrix[i][k] * factor;
        if (k < j) {
          lowerTriangularMatrix[j][k] = factor;
        }
      }
    }
  }
  
  /** 
   * We know that ax = b, after decompose a into l and u matrices, we get lux = b
   * Supose that ux = y and ly = b.
   */

  for (let i = 0; i < dimension; i++) {
    lowerTriangularMatrix[i].push(system[i][dimension]);
  }

  const y = lowerTriangularSystemSolver(lowerTriangularMatrix);

  for (let i = 0; i < dimension; i++) {
    upperTriangularMatrix[i][dimension] = y[i];
  }

  return upperTriangularSystemSolver(upperTriangularMatrix);
}

/**************************************************************************************************
 *                                          Jacobi method                                         *
 **************************************************************************************************/

/**
 * This function used to solve a linear system of equations using Jacobi method
 * @param {Array} system
 * @returns {Array|Boolean} 
 */

function jacobiMethod(system) {
  if (isApplicableJacobiMethod(system) == false) return false;

  const numberOfRows = system.length;
  const numberOfCols = system[0].length;
  const resultVector = new Array(numberOfRows);
}

/**
 * Method to check if the system is a valid jacobi system
 * @param {Array} system 
 * @returns {Boolean}
 */

function isApplicableJacobiMethod(system) {
  const numberOfRows = system.length;

  for (let i = 0; i < numberOfRows; i++) {
    if (system[i][i] == 0) {
      return false;
    }
  }

  return true;
}

/**
 * Function to check wethere the system matrix is strictlly dominante
 * @param {Array}
 * @returns {Boolean}
 */

function isStrictlyDominante(system) {
  const dimension = system.length;

  for (let i = 0; i < dimension; i++) {
    const compared = Math.abs(system[i][i]);
    var sum1 = system[i][0], isGreater1 = false;
    var sum2 = system[0][i], isGreater2 = false;

    for (let j = 1; j < dimension && !isGreater1 && !isGreater2; j++) {
      sum1 += i != j ? Math.abs(system[i][j]):0;
      sum2 += i != j ? Math.abs(system[j][i]):0;
      isGreater1 = sum1 >= compared;
      isGreater2 = sum2 >= compared;
    }

    if (isGreater1 && isGreater2) {
      return false;
    }
  }

  return true;
}
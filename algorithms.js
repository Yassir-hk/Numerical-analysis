/**************************************************************************************************
 *                                              Handlers                                          *
 **************************************************************************************************/

const regexOpt = /^[+-]$/;
const regexNum = /^[\d.]+$/;
const regexVar = /^[a-zA-Z]+$/;

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
 * @param {Array} equation
 * @returns {Boolean}
 */

function isValidEquation(equation) {
  let lastChar = false;
  let countEquales = 0;

  for (let char of equation) {
    const isValidNum = regexNum.test(char) && (regexOpt.test(lastChar) || regexNum.test(lastChar) || lastChar == false);
    const isValidOpt = regexOpt.test(char) && (regexVar.test(lastChar) || lastChar == false || lastChar == '=');
    const isValidVar = regexVar.test(char) && regexNum.test(lastChar);

    if (!isValidNum && !isValidOpt && !isValidVar) {
      return false;
    }

    countEquales += (char == '=');
    lastChar = char;
  }

  return countEquales == 1;
}

/**
 * Function to verity the equations input.
 * @returns {Array | boolean}
 */

function fetchSystem() {
  const equations = Array.from(document.getElementsByClassName("equationInput"));
  const variablesOfSystemSet = new Set();
  const system = [];

  for (let equation of equations) {
    const equationString = equation.value.replace(/\s+/g, '');

    // Check if the equation is valid.
    // if (!isValidEquation(equationString)) {
    //   return false;
    // }

    const systemEquation = {};
    let coefficient = '', operator = '';

    for (let char of equationString) {
      if (regexNum.test(char)) {
        coefficient += char;
      } else if (regexVar.test(char)) {
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
      } else if (regexOpt.test(char)) {
        operator = char;
      }
    }
    if (coefficient.length) {
      systemEquation['='] = parseFloat((operator == '=' ? '' : operator) + coefficient);
    }

    system.push(systemEquation);
  }

  // Check if the number of variables equales the number of equations.
  if (variablesOfSystemSet.size != equations.length) {
    return false;
  }

  // Add coefficients with value 0.
  const variablesOfSystemArray = [...variablesOfSystemSet].sort();
  const matrixOfSystem = [];

  for (let systemEquation of system) {
    const equationCoefficients = [];

    for (let variable of variablesOfSystemArray) {
      if (variable in systemEquation) {
        equationCoefficients.push(systemEquation[variable]);
      } else {
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
 * @returns {Boolean|void}
 */

function compute() {
  const system = fetchSystem();
  const computationMethod = parseInt(document.getElementById("computationMethod").value);

  // Check if the system is valid.
  if (system == false) {
    alert("This system cannot be processed");
    return false;
  }

  const resultVectorContainer = document.querySelector(".container-2");
  resultVectorContainer.style.visibility = "visible";

  // Use the selected method to compute the result.
  switch (computationMethod) {
    case 1:
      writeResult(gaussElemination(system[0]), system[1]);
      break;
    case 2:
      writeResult(luDecomposition(system[0]), system[1]);
      break;
    case 3:
      writeResult(jacobiMethod(system[0]), system[1]);
      break;
    case 4:
      writeResult(gaussSeidelMethod(system[0]), system[1]);
      break;
  }
}

/**
 * Function to write the solution of the system in the resultVector container.
 * @param {Array | Boolean} resultVector
 * @param {Array} variablesOfSystem
 * @returns {Boolean}
 */

function writeResult(resultVector, variablesOfSystem) {
  const resultVectorContainer = document.getElementById("resultContainer");

  if (resultVector == false) {
    resultVectorContainer.innerHTML = "This system has infinity of solutions or no solution.";
    return false;
  }

  let resultVectorContent = '';
  let currentVariableIndx = 0;

  for (let value of resultVector) {
    resultVectorContent += variablesOfSystem[currentVariableIndx++] + " = " + value + "<br><br>";
  }

  resultVectorContainer.innerHTML = resultVectorContent;
}

/**************************************************************************************************
 *                                      Matrice functions                                         *
 **************************************************************************************************/

/**
 * Function to multiply two matrices.
 * @param {Array} firstMatrix
 * @param {Array} secondMatrix
 * @returns {Array | Boolean}
 */

function multiplyTwoMatrices(firstMatrix, secondMatrix) {
  if (firstMatrix[0].length != secondMatrix.length) {
    return false;
  }

  const numberOfRows = firstMatrix.length;
  const numberOfCols = secondMatrix[0].length;
  const resultMatrix = new Array(numberOfRows);

  for (let i = 0; i < numberOfRows; i++) {
    resultMatrix[i] = new Array(numberOfCols);

    for (let j = 0; j < numberOfCols; j++) {
      for (let k = 0; k < numberOfCols; k++) {
        resultMatrix[i][j] += (firstMatrix[i][k] * secondMatrix[k][j]);
      }
    }
  }

  return resultMatrix;
}

/**
 * Function to calculate the addition of two matrices.
 * @param {Array} firstMatrix
 * @param {Array} secondMatrix
 * @returns {Array | Boolean}
 */

function addTwoMatrices(firstMatrix, secondMatrix) {
  if (firstMatrix.length != secondMatrix.length || firstMatrix[0].length != secondMatrix.length) {
    return false;
  }

  const numberOfRows = firstMatrix.length;
  const numberOfCols = firstMatrix[0].length;
  const resultMatrix = new Array(numberOfRows);

  for (let i = 0; i < numberOfRows; i++) {
    resultMatrix[i] = new Array(numberOfCols);

    for (let j = 0; j < numberOfCols; j++) {
      resultMatrix[i][j] = firstMatrix[i][j] + secondMatrix[i][j];
    }
  }

  return resultMatrix;
}

/**
 * Function to check wethere the system matrix is strictlly diagonally dominante 
 * @param {Array}
 * @returns {Boolean}
 */

function isStrictlyDiagonallyDominant(system) {
  const dimension = system.length;

  for (let i = 0; i < dimension; i++) {
    const diagonalElement = Math.abs(system[i][i]);

    if (diagonalElement == 0) {
      return false;
    }

    for (let j = 0, rowSum = 0; j < dimension; j++) {
      rowSum += (i != j ? Math.abs(system[i][j]):0);
      if (rowSum >= diagonalElement) {
        return false;
      }
    }
  }

  return true;
}

/**************************************************************************************************
 *                                      Gauss elemination                                         *
 **************************************************************************************************/

/**
 * Function is used to solve a system of linear equations using the Gaussian elimination method.  
 * @param {Array} system 
 * @returns {Array | Boolean}
 */

function gaussElemination(system) {
  const numberOfRows = system.length;
  const numberOfCols = system[0].length;

  // Triangularization of the system.
  for (let i = 0; i < numberOfRows; i++) {
    let pivotIndex = i;

    while (pivotIndex < numberOfRows && system[pivotIndex][i] == 0) {
      pivotIndex++;
    }

    if (pivotIndex < numberOfRows && pivotIndex != i) {
      [system[pivotIndex], system[i]] = [system[i], system[pivotIndex]];
    } else if (pivotIndex == numberOfRows) {
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
 * Function that solve the upper triangular system.
 * @param {Array} system
 * @returns {Array | Boolean}
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
 * Function that solve the lower triangular system.
 * @param {Array} system
 * @returns {Array | Boolean}
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
 *                                        LU decomposition                                        *
 **************************************************************************************************/

/**
 * Function to get the lower and the upper matrices.
 * @param {Array} system
 * @returns {Array | Boolean}
 */

function luDecomposition(system) {
  const dimension = system.length;
  const lowerTriangularMatrix = [];
  const upperTriangularMatrix = Array.from(system);

  // Initializing the lower matrix.  
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
    } else if (pivotIndex == dimension) {
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
   * We know that Ax = b, after decompose A into L and U matrices, we get LUx = b.
   * Supose that Ux = y and Ly = b.
   * Step 1: solve Ly = b.
   * Step 2: solve Ux = y.
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
 * This function used to solve a linear system of equations using Jacobi method.
 * @param {Array} system
 * @returns {Array | Boolean} 
 */

function jacobiMethod(system) {
  if (!isStrictlyDiagonallyDominant(system)) {
    return false;
  }

  const numberOfRows = system.length;
  const numberOfCols = system[0].length;
  const resultVector = new Array(numberOfRows).fill(0);

  let iterationsLimit = 1e3;
  let estimationError = 1e3;

  while (estimationError >= 1 && iterationsLimit--) {
    for (let i = 0; i < numberOfRows; i++) {
      let sum = 0;

      for (let j = 0; j < numberOfCols - 1; j++) {
        if (i == j) {
          continue;
        }
        sum += system[i][j] * resultVector[j];
      }

      const prevResultVector = resultVector[i];
      resultVector[i] = (system[i][numberOfCols - 1] - sum) / system[i][i];
      estimationError = Math.abs(resultVector[i] - prevResultVector);
    }
  }

  return resultVector;
}

/**************************************************************************************************
 *                                          Gauss-Seidel method                                   *
 **************************************************************************************************/

/**
 * This function solve a linear system of equations using the Gauss-Seidel method which is an extension of the Jacobi method.
 * @param {Array} system
 * @returns {Array | Boolean}
 */

function gaussSeidelMethod(system) {
  if (!isStrictlyDiagonallyDominant(system)) {
    return false;
  }

  const numberOfRows = system.length;
  const numberOfCols = system[0].length;
  const resultVector = new Array(numberOfRows).fill(0);

  let iterationsLimit = 1e3;
  let estimationError = 1e3;

  while (estimationError >= 1 && iterationsLimit--) {
    for (let i = 0; i < numberOfRows; i++) {
      let sum = 0;

      for (let j = 0; j < numberOfRows; j++) {
        if (i == j) {
          continue;
        }
        const product = resultVector[j] * system[i][j];
        sum += product * (j < i ? -1:1);
      }

      const prevResultVector = resultVector[i];
      resultVector[i] = (system[i][numberOfCols - 1] - sum) / system[i][i];
      estimationError = Math.abs(resultVector[i] - prevResultVector);
    }
  }

  return resultVector;
}

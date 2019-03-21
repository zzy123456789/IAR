package IAR;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.Arrays;

import common.*;

/**
 * 
 * @author ziyin
 *
 */
public class IAR {
	/**
	 * The item number
	 */
	int itemNumber;

	/**
	 * The item information
	 */
	ItemInformation[][] itemMatrix;

	/**
	 * The item degree
	 */
	double[] itemDegree;

	/**
	 * The information for leaving one out
	 */
	LeaveOneInformation leaveOneInfo;

	public static double[] M;// the quality of item
	
	public static double MAE;
	public static double RMSE;
	
	public static int max;

	/**
	 * 
	 * @param paraFileName
	 * @throws Exception
	 */
	public IAR(String paraFileName, int paraUserNum, int paraItemNum) throws Exception {
		File tempFile = null;

		String tempString = null;
		double tempRating = 0;
		int tempUserIndex = 0;
		int tempItemIndex = 0;
		int tempTime = 0;
		String[] tempStrArray = null;
		M = new double[paraItemNum];

		// Compute values of arrays
		tempFile = new File(paraFileName);
		if (!tempFile.exists()) {
			System.out.println("File is not exist!");
			return;
		} // Of if

		RandomAccessFile tempRanFile = new RandomAccessFile(tempFile, "r");
		// 锟斤拷锟侥硷拷锟斤拷锟斤拷始位锟斤拷
		int tempBeginIndex = 0;
		// 锟斤拷锟斤拷锟侥硷拷锟侥匡拷始位锟斤拷锟狡碉拷beginIndex位锟矫★拷
		tempRanFile.seek(tempBeginIndex);

		// Step 1. count the item degree
		int[] tempItemDegree = new int[paraItemNum];
		while ((tempString = tempRanFile.readLine()) != null) {
			tempStrArray = tempString.split(",");
			tempItemIndex = Integer.parseInt(tempStrArray[1]);

			tempItemDegree[tempItemIndex]++;
		} // Of while
		
		max = 0;

		for (int i = 0; i < paraItemNum; i++) {
			if(tempItemDegree[i] > max) {
				max = tempItemDegree[i];
			}
			// System.out.println(tempItemDegree[i]);
			// System.out.println(paraItemNum);

			// //M1
			 if(tempItemDegree[i] == 0) {
			 continue;
			 }
			 M[i] = (double) (Math.log10((double) tempItemDegree[i]));
//			 System.out.println(M[i]);

			// M2

//			 if (tempItemDegree[i] == 0) {
//			 continue;
//			 }
//			 if (tempItemDegree[i] == paraUserNum) {
////			 System.out.println("the log"
////			 + 2 * 1 / (double) Math.abs((Math.log10((double) (paraUserNum - 1) / (double)
////			 paraUserNum))));
//			 M[i] = 2 * 1 / (double) Math.abs((Math.log10((double) (paraUserNum - 1) /
//			 (double) paraUserNum)));
//			 } else {
//			 M[i] = 1.0
//			 / (double) Math.abs((double) (Math.log10((double) tempItemDegree[i] /
//			 (double) paraUserNum)));
//			 }
//			 System.out.println(M[i]);
			//
			// //M3
//			if (tempItemDegree[i] == 0) {
//				continue;
//			}
//			if (tempItemDegree[i] == paraUserNum) {
////				System.out.println(2*1/Math.log10(paraUserNum / (double)(paraUserNum - 1)));// (Math.log(paraUserNum /(paraUserNum
//																				// -1)))/(double)Math.log(10)
//				M[i] = (double)tempItemDegree[i]*2*1/Math.log10(paraUserNum / (double)(paraUserNum - 1));
//			}else {
//			M[i] = (double)tempItemDegree[i]*1.0/ (double) Math.abs((double) (Math.log10((double) tempItemDegree[i] /
//					 (double) paraUserNum)));
//			}
//			System.out.println(M[i]);
		}

		// Step 2. Construct the item information
		// Step 2.1 Allocate the memory for the item information
		ItemInformation[][] tempItemMatrix = new ItemInformation[paraItemNum][];
		for (int i = 0; i < paraItemNum; i++) {
			tempItemMatrix[i] = new ItemInformation[tempItemDegree[i]];
		} // Of for i

		for (int i = 0; i < tempItemDegree.length; i++) {
			tempItemDegree[i] = 0;
		} // Of for i

		// Step 2.2 Reread the file and construct the item information
		tempRanFile.seek(tempBeginIndex);
		while ((tempString = tempRanFile.readLine()) != null) {
			// System.out.println(tempString);
			tempStrArray = tempString.split(",");
			tempUserIndex = Integer.parseInt(tempStrArray[0]);
			tempItemIndex = Integer.parseInt(tempStrArray[1]);
			tempRating = Double.parseDouble(tempStrArray[2]);
			// tempTime = Integer.parseInt(tempStrArray[3]);

			tempItemMatrix[tempItemIndex][tempItemDegree[tempItemIndex]] = new ItemInformation();
			tempItemMatrix[tempItemIndex][tempItemDegree[tempItemIndex]].userIndex = tempUserIndex;
			tempItemMatrix[tempItemIndex][tempItemDegree[tempItemIndex]].itemIndex = tempItemIndex;
			tempItemMatrix[tempItemIndex][tempItemDegree[tempItemIndex]].rating = tempRating;
			// tempItemMatrix[tempItemIndex][tempItemDegree[tempItemIndex]].time = tempTime;
			tempItemDegree[tempItemIndex]++;
			// System.out.println(tempTime);
		} // Of while
		tempRanFile.close();

		// Step 2.3 Count the item number
		for (int i = 0; i < tempItemDegree.length; i++) {
			if (tempItemDegree[i] > 1e-6) {
				itemNumber++;
			} // Of if
		} // Of for i

		// Step 2.4 Compress
		buildItemMatrix(tempItemMatrix, tempItemDegree);
	}// Of the first constructor

	/**
	 * Build the item-based information matrix
	 */
	void buildItemMatrix(ItemInformation[][] paraItemInfo, int[] paraItemDegree) {
		itemMatrix = new ItemInformation[itemNumber][];
		itemDegree = new double[itemNumber];
		int tempCompressIndex = 0;
		for (int i = 0; i < paraItemDegree.length; i++) {
			if (paraItemDegree[i] < 1e-6) {
				continue;
			} // Of if
			itemMatrix[tempCompressIndex] = new ItemInformation[paraItemDegree[i]];
			for (int j = 0; j < paraItemDegree[i]; j++) {
				itemMatrix[tempCompressIndex][j] = new ItemInformation();
				itemMatrix[tempCompressIndex][j].userIndex = paraItemInfo[i][j].userIndex;
				itemMatrix[tempCompressIndex][j].itemIndex = paraItemInfo[i][j].itemIndex;
				itemMatrix[tempCompressIndex][j].rating = paraItemInfo[i][j].rating;
				// itemMatrix[tempCompressIndex][j].time = paraItemInfo[i][j].time;
			} // Of for j
			itemDegree[tempCompressIndex] = paraItemDegree[i];
			tempCompressIndex++;
		} // Of for i

		// printMatrix(itemMatrix);

	}// Of buildItemMatrix

	/**
	 * Leave one out
	 */
	void leaveOneOut(int paraX, int paraY, double paraRating) {
		leaveOneInfo = new LeaveOneInformation();
		leaveOneInfo.x = paraX;
		leaveOneInfo.y = paraY;
		leaveOneInfo.leaveRating = paraRating;
		leaveOneInfo.userIndex = itemMatrix[paraX][paraY].userIndex;
		leaveOneInfo.itemIndex = itemMatrix[paraX][paraY].itemIndex;

		itemMatrix[paraX][paraY].rating = 0;
		itemDegree[paraX]--;
	}// of leaveOneOut

	/**
	 * Leave one out restore
	 */
	void leaveOneOutRestore() {
		itemMatrix[leaveOneInfo.x][leaveOneInfo.y].rating = leaveOneInfo.leaveRating;

		itemDegree[leaveOneInfo.x]++;
	}// Of leaveOneOutRestore

	/**
	 * prediction
	 */
	double prediction(int paraX, int paraY, int paraK) {
		// step 1. leave one out
		double tempRating = itemMatrix[paraX][paraY].rating;
		// System.out.println(tempRating);
		leaveOneOut(paraX, paraY, tempRating);
		// System.out.println("leaveOneOut: " + paraX + " : " + paraY + " = " +
		// tempRating);

		// Step 2. compute the similarity
		// Step 2.1 initialize the distances of neighbors
		double[] tempNeighborDistances = new double[paraK + 2];
		double[] tempNeighborRating = new double[paraK + 2];
		for (int i = 0; i < tempNeighborRating.length; i++) {
			tempNeighborRating[i] = -1;
		} // Of for
		for (int i = 1; i <= paraK + 1; i++) {
			tempNeighborDistances[i] = -Double.MAX_VALUE;
		} // Of for i
		tempNeighborDistances[0] = Double.MAX_VALUE;// A sign
		double tempDistance;
		// System.out.println("Current neiborghood array:" +
		// Arrays.toString(tempNeighborDistances));

		// Step 2.2 Find the neighbors
		for (int i = 0; i < itemMatrix.length; i++) {
			// step 2.2.1 exclude myself
			if (i == paraX) {
				continue;
			} // Of if
				// System.out.println("i: " + i);
				// Step 2.2.2 search the user index
			for (int j = 0; j < itemMatrix[i].length; j++) {
				// System.out.println("test: " + leaveOneInfo.userIndex + " : " +
				// itemMatrix[i][j].userIndex);
				if (leaveOneInfo.userIndex == itemMatrix[i][j].userIndex) {
					// System.out.println(i + " : " + paraX);
					tempDistance = gravityDistance(i, paraX, leaveOneInfo.userIndex);
					// tempDistance = euclideanDistance(i, paraX);
					// tempDistance = manhattanDistance(i, paraX);
//					 System.out.println("tempDistance = " + tempDistance);

					if (tempDistance == 0) {
						continue;
					}

					// Insert(锟斤拷小锟斤拷锟斤拷锟斤拷锟斤拷)
					for (int k = paraK; k >= 0; k--) {
//						 System.out.println("tempDistance: k : neigborDistance " +
//						 tempDistance +" ; " + k + " ; " + tempNeighborDistances[k]);
//						 System.out.println("k = " + k);
						if (tempDistance > tempNeighborDistances[k]) {
							if (k == 0) {
//								 System.out.println("Distance: " + tempDistance);
							} // Of if
							tempNeighborDistances[k] = tempNeighborDistances[k - 1];
							tempNeighborRating[k] = tempNeighborRating[k - 1];
//							 System.out.println("tempNeighborDistances1[k] = " +
//							 tempNeighborDistances[k]);
						} else {// if(tempDistance <= tempNeighborDistances[k])//引力
							tempNeighborDistances[k +1] = tempDistance;
							tempNeighborRating[k +1] = itemMatrix[i][j].rating; //
//							 System.out.println("tempNeighborDistances2[k] = " +
//							 tempNeighborDistances[k]);
							break;
						} // Of if
					} // Of for j
//						 System.out.println("Current neiborghood array:" +
//						 Arrays.toString(tempNeighborDistances));
				} // Of if
			} // Of for j
		} // Of for i

		// Step 2.3 Compute the prediction
		int tempCount = 0;
		double tempPrediction = 0;
		tempNeighborRating[paraK + 1] = -1;
		for (int i = 1; i < paraK + 1; i++) {
			if (tempNeighborRating[i] != -1) {
				tempPrediction += tempNeighborRating[i];
				tempCount++;
			} // Of if
		} // Of for i
			// Step 3. leave one out restore
		leaveOneOutRestore();

		if (tempCount != 0) {
			tempPrediction /= tempCount;
		} // Of if

		return tempPrediction;
	}// Of prediction

	/**
	 * Compute MAE
	 */
	void computeMAEandRMSE(int paraK) {
		double tempMAE = 0;
		double tempRMSE = 0;
		double tempPrediction = 0;
		double tempCount = 0;
		// System.out.println("MAE");
		for (int i = 0; i < itemMatrix.length; i++) {
			for (int j = 0; j < itemMatrix[i].length; j++) {

				// System.out.println("the i : j is " + i + " : " + j);
				// if(itemMatrix[i][j].rating == 0) {
				// continue;
				// }
				tempPrediction = prediction(i, j, paraK);

//				 System.out.println("the prediction is " + tempPrediction);

				tempMAE += Math.abs(itemMatrix[i][j].rating - tempPrediction);
				tempRMSE += Math.pow((itemMatrix[i][j].rating - tempPrediction), 2);

				tempCount++;
			} // Of for j
		} // Of for i

		MAE =  tempMAE / tempCount;
		RMSE = Math.sqrt(tempMAE / tempCount);
	}// Of computeMAE

	/**
	 * Gravity distance
	 */
	double gravityDistance(int paraIndex1, int paraIndex2, int paraUserIndex) {
		double tempDistance = 0;
		double result = 0;
		// int tempCount = 0;
		// double tempRating1 = 0;
		// double tempRating2 = 0;
		// double tempRating = 0;

		// for (int i = 0; i < itemMatrix[paraIndex1].length; i++) {
		// for (int j = 0; j < itemMatrix[paraIndex2].length; j++) {
		// if (itemMatrix[paraIndex1][i].userIndex ==
		// itemMatrix[paraIndex2][j].userIndex) {
		// tempRating1 = itemMatrix[paraIndex1][i].rating;
		// tempRating2 = itemMatrix[paraIndex2][j].rating;
		// // if (tempRating1 != 0 && tempRating2 != 0) {
		// tempRating += Math.pow((tempRating1 - tempRating2), 2);
		// tempCount++;
		// } // Of if
		// } // Of for j
		// } // Of for i
		//
		// if (tempCount == 0) {
		// return 0;
		// }

		// tempDistance = tempRating / tempCount;
//		tempDistance = euclideanDistance(paraIndex1, paraIndex2, paraUserIndex);
		 tempDistance = manhattanDistance(paraIndex1, paraIndex2, paraUserIndex);
//		 System.out.println("the distance is " + tempDistance);

		if ((M[paraIndex1] * M[paraIndex2]) == 0) {
			return 0;
		}

		if (tempDistance == 0) {
			return Integer.MAX_VALUE;
		}
		
//		System.out.println(M[paraIndex1] + " : " + M[paraIndex2]);

		result = (double) (M[paraIndex1] * M[paraIndex2]) / Math.pow(tempDistance, 2);

		return result;
	}// Of gravityDistance

	/**
	 * Manhattan distance
	 */
	double manhattanDistance(int paraIndex1, int paraIndex2, int paraUserIndex) {
		double tempDistance = 0;
		double result = 0;
		int tempCount = 0;
		double tempRating1 = 0;
		double tempRating2 = 0;
		double tempRating = 0;

		for (int i = 0; i < itemMatrix[paraIndex1].length; i++) {
			for (int j = 0; j < itemMatrix[paraIndex2].length; j++) {
				if (itemMatrix[paraIndex1][i].userIndex == itemMatrix[paraIndex2][j].userIndex
						&& itemMatrix[paraIndex1][i].userIndex != paraUserIndex) {
					tempRating1 = itemMatrix[paraIndex1][i].rating;
					tempRating2 = itemMatrix[paraIndex2][j].rating;
					// if (tempRating1 != 0 && tempRating2 != 0) {
					tempRating += Math.abs(tempRating1 - tempRating2);
//					 System.out.println(Math.abs(tempRating1 - tempRating2));
					tempCount++;
				} // Of if
			} // Of for j
		} // Of for i

		if (tempCount == 0) {
			return 0;
		}
		return tempRating;
	}// Of gravityDistance

	/**
	 * Euclidean distance
	 */
	double euclideanDistance(int paraIndex1, int paraIndex2, int paraUserIndex) {
		double tempDistance = 0;
		double result = 0;
		int tempCount = 0;
		double tempRating1 = 0;
		double tempRating2 = 0;
		double tempRating = 0;

		for (int i = 0; i < itemMatrix[paraIndex1].length; i++) {
			for (int j = 0; j < itemMatrix[paraIndex2].length; j++) {
				if (itemMatrix[paraIndex1][i].userIndex == itemMatrix[paraIndex2][j].userIndex
						&& itemMatrix[paraIndex1][i].userIndex != paraUserIndex) {
					tempRating1 = itemMatrix[paraIndex1][i].rating;
					tempRating2 = itemMatrix[paraIndex2][j].rating;
					// if (tempRating1 != 0 && tempRating2 != 0) {
					tempRating += Math.pow((tempRating1 - tempRating2), 2);
//					 System.out.println(Math.pow((tempRating1 - tempRating2), 2));
					tempCount++;
				} // Of if
			} // Of for j
		} // Of for i

		if (tempCount == 0) {
			return 0;
		}

		result = Math.sqrt(tempRating);

		return result;
	}// Of gravityDistance

	/**
	 * 
	 * @param paraMatrix
	 */
	public void printMatrix(ItemInformation[][] paraMatrix) {
		for (int i = 0; i < paraMatrix.length; i++) {
			if (paraMatrix[i] == null) {
				continue;
			} // Of if
			for (int j = 0; j < paraMatrix[i].length; j++) {
				System.out.print(paraMatrix[i][j].userIndex + "," + paraMatrix[i][j].itemIndex + ","
						+ paraMatrix[i][j].rating + ";");
			} // Of for j
			System.out.println();
		} // Of for i
	}// Of printMatrix

	/**
	 * 
	 * @param paraMatrix
	 */
	public void printMatrix(boolean[][] paraMatrix) {
		for (int i = 0; i < paraMatrix.length; i++) {
			if (paraMatrix[i] == null) {
				continue;
			} // Of if
			for (int j = 0; j < paraMatrix[i].length; j++) {
				System.out.print(paraMatrix[i][j] + ";");
			} // Of for j
			System.out.println();
		} // Of for i
	}// Of printMatrix

	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO 锟皆讹拷锟斤拷锟缴的凤拷锟斤拷锟斤拷锟�
		try {
			IAR tempTime = new IAR("data/movielens943u1682m.txt", 943, 1682);// "data/test.txt", 4, 5  "data/test1.txt", 5, 3
															// "data/movielens943u1682m.txt", 943, 1682

			// tempTime.printMatrix(tempTime.itemMatrix);
			// System.out.println("The g is" + Arrays.toString(tempTime.itemDegree));
			// System.out.println("The M is" + Arrays.toString(tempTime.M));
//			 for (int i = 10; i <=100; i+=10) {
//
////			int i = 10;
//			tempTime.computeMAEandRMSE(i);
//			System.out.println("MAE: " + tempTime.MAE + " RMSE: " + tempTime.RMSE);
//
//			 }
			System.out.println(tempTime.max);

		} catch (Exception ee) {
			ee.printStackTrace();
		} // Of try
	}// Of main

}// Of class TimeSeriesRecommendation

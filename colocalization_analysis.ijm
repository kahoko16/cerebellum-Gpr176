macro "Colocalization Analysis (Receptor vs Lysosome)" {

    // ============================================================
    // 受容体(緑)チャンネルとリソソーム(赤)チャンネルの
    // colocalization (Pearson's R, Manders' M1/M2) を
    // 細胞ROIごとに計算するマクロ。
    //
    // intensity_analysis.ijm と同じROI(.zip)をそのまま使い回せる
    // （ROIは形状情報なのでチャンネル・解析内容に依存しない）。
    //
    // 前提:
    // - 緑チャンネル・赤チャンネルは同じファイル名(basename)で
    //   別々のフォルダに16bit tifとして保存されている
    // - ROI(.zip)も同じbasenameで保存されている
    // - ROI Managerの最後のROIは背景ROIとして解析から除外する
    //   (intensity_analysis.ijm と同じ約束)
    //
    // 注意: 実際のファイルでの動作確認はしていません。
    // 必ず1〜2ファイルで試してから一括実行してください。
    // ============================================================

    // 各チャンネルの背景値(Fbck)。background測定で求めた値を入れてください。
    // Manders' M1/M2 の閾値として使用します。
    greenBackground = 0;  // 受容体(緑)チャンネルの背景値
    redBackground = 0;    // リソソーム(赤)チャンネルの背景値

    greenDir = getDirectory("受容体(緑)チャンネルの16bit tifが入っているフォルダを選択してください");
    if (greenDir == "") exit("フォルダが選択されませんでした。");

    redDir = getDirectory("リソソーム(赤)チャンネルの16bit tifが入っているフォルダを選択してください");
    if (redDir == "") exit("フォルダが選択されませんでした。");

    roiDir = getDirectory("ROI(.zip)が入っているフォルダを選択してください");
    if (roiDir == "") exit("フォルダが選択されませんでした。");

    outputDir = getDirectory("結果を保存するフォルダを選択してください");
    if (outputDir == "") exit("フォルダが選択されませんでした。");

    list = getFileList(greenDir);
    setBatchMode(true);
    run("Clear Results");

    for (f = 0; f < list.length; f++) {
        fileName = list[f];
        lowerName = toLowerCase(fileName);
        if (!endsWith(lowerName, ".tif") && !endsWith(lowerName, ".tiff")) continue;

        dotIndex = lastIndexOf(fileName, ".");
        basename = substring(fileName, 0, dotIndex);

        redPath = redDir + fileName;
        if (!File.exists(redPath)) {
            print("スキップ (対応する赤チャンネル画像が見つかりません): " + fileName);
            continue;
        }

        roiPath = roiDir + basename + ".zip";
        if (!File.exists(roiPath)) {
            print("スキップ (対応するROIが見つかりません): " + fileName);
            continue;
        }

        open(greenDir + fileName);
        greenID = getImageID();

        open(redPath);
        redID = getImageID();

        roiManager("reset");
        roiManager("open", roiPath);
        nTotal = roiManager("count");
        if (nTotal < 1) {
            selectImage(greenID); close();
            selectImage(redID); close();
            continue;
        }

        // 背景ROI(最後のROI)は解析から除外する（intensity_analysis.ijmと同じ約束）
        bgIndex = nTotal - 1;

        for (i = 0; i < nTotal; i++) {
            if (i == bgIndex) continue;

            selectImage(greenID);
            roiManager("select", i);
            Roi.getContainedPoints(xpts, ypts);
            n = xpts.length;
            if (n < 1) continue;

            greenVals = newArray(n);
            redVals = newArray(n);

            for (k = 0; k < n; k++) {
                selectImage(greenID);
                greenVals[k] = getPixel(xpts[k], ypts[k]);
                selectImage(redID);
                redVals[k] = getPixel(xpts[k], ypts[k]);
            }

            pearsonR = computePearson(greenVals, redVals, n);

            // Manders' M1/M2 (各チャンネルの背景値を閾値として使用)
            greenSum = 0;
            redSum = 0;
            greenAboveRedThresh = 0;
            redAboveGreenThresh = 0;
            for (k = 0; k < n; k++) {
                g = greenVals[k];
                r = redVals[k];
                greenSum += g;
                redSum += r;
                if (r > redBackground) greenAboveRedThresh += g;
                if (g > greenBackground) redAboveGreenThresh += r;
            }
            m1 = 0;
            if (greenSum > 0) m1 = greenAboveRedThresh / greenSum;
            m2 = 0;
            if (redSum > 0) m2 = redAboveGreenThresh / redSum;

            row = nResults;
            setResult("Image_Name", row, basename);
            setResult("Cell_ID", row, i + 1);
            setResult("N_Pixels", row, n);
            setResult("Pearson_R", row, pearsonR);
            setResult("Manders_M1_GreenInRed", row, m1);
            setResult("Manders_M2_RedInGreen", row, m2);
        }

        selectImage(greenID); close();
        selectImage(redID); close();
    }

    updateResults();
    saveAs("Results", outputDir + "colocalization_results.csv");
    setBatchMode(false);
    showMessage("Colocalization解析が完了しました！");
}

function computePearson(a, b, n) {
    sumA = 0;
    sumB = 0;
    for (k = 0; k < n; k++) {
        sumA += a[k];
        sumB += b[k];
    }
    meanA = sumA / n;
    meanB = sumB / n;

    num = 0;
    denomA = 0;
    denomB = 0;
    for (k = 0; k < n; k++) {
        da = a[k] - meanA;
        db = b[k] - meanB;
        num += da * db;
        denomA += da * da;
        denomB += db * db;
    }
    if (denomA == 0 || denomB == 0) return 0;
    return num / sqrt(denomA * denomB);
}

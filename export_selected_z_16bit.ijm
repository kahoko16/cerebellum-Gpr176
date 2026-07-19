macro "Export Selected Z as 16-bit from Olympus files" {

    // ============================================================
    // 「良いZ」として選んで残してあるtif/zip(8bit)のファイル名から
    // チャンネル・Z番号を読み取り、対応する元のOlympusファイル
    // (.oib / .oir / .oif) からそのチャンネル・Zスライスだけを16bitのまま
    // 読み込み直して、同じファイル名でtif保存するマクロ。
    //
    // 出力ファイル名は入力と同じ basename になるので、既存の
    // ROI (.zip) をそのまま適用できる（ROIは座標情報なのでbit深度に
    // 依存しない）。実行するたびに指定フォルダの下にタイムスタンプ付き
    // サブフォルダ(run_YYYYMMDD_HHMMSS)を新規作成してそこに書き出すので、
    // 過去の出力を上書き・削除することはない。
    //
    // 前提: ファイル名が "○○○_C001Z004.tif" のように
    // 末尾に "_C<チャンネル3桁>Z<Z番号3桁>" が付いている形式。
    // 元のOlympusファイルは "○○○.oib" "○○○.oir" "○○○.oif" のいずれか。
    //
    // 注意: 実際のファイルでの動作確認はしていません。
    // 必ず1〜2ファイルで試してから一括実行してください。
    // ============================================================

    // 定量に使うチャンネル番号(Olympusファイル内のチャンネル通し番号)。
    // ファイル名の "_C001" 等の数字とは別物の場合があるので要確認。
    targetChannel = 1;

    // tif作成時に元のOlympusファイル名には無い波長サフィックスを
    // 追加している場合、ここに指定すると元ファイル名から取り除く
    // (例: "260703-WT-DMSO-before-1206_02-488" -> "260703-WT-DMSO-before-1206_02")
    // ファイルによって "-488" と "=488" のように表記がゆれている場合は、
    // 両方とも配列に入れておけば両方試す。不要なら空配列 newArray() のままにする
    wavelengthSuffixCandidates = newArray("-488", "=488");

    selectedDir = getDirectory("選んだZのtif/zipが入っているフォルダを選択してください");
    if (selectedDir == "") exit("フォルダが選択されませんでした。");

    olympusDir = getDirectory("元のOlympusファイル(.oib/.oir/.oif)が入っているフォルダを選択してください");
    if (olympusDir == "") exit("フォルダが選択されませんでした。");

    outputDirBase = getDirectory("16bit tifの保存先フォルダを選択してください");
    if (outputDirBase == "") exit("フォルダが選択されませんでした。");

    // 実行するたびにタイムスタンプ付きのサブフォルダに書き出すので、
    // 既存の出力を上書き/削除することはない（basenameは変えないので
    // ROI(.zip)との対応関係はそのフォルダ内でそのまま使える）
    getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
    runFolderName = "run_" + year + IJ.pad(month + 1, 2) + IJ.pad(dayOfMonth, 2) +
                     "_" + IJ.pad(hour, 2) + IJ.pad(minute, 2) + IJ.pad(second, 2);
    outputDir = outputDirBase + runFolderName + "/";
    File.makeDirectory(outputDir);
    print("出力先フォルダ: " + outputDir);

    list = getFileList(selectedDir);
    nDone = 0;
    nSkipped = 0;

    for (i = 0; i < list.length; i++) {
        fileName = list[i];
        lowerName = toLowerCase(fileName);
        if (!endsWith(lowerName, ".tif") && !endsWith(lowerName, ".tiff")) continue;

        dotIndex = lastIndexOf(fileName, ".");
        basename = substring(fileName, 0, dotIndex);

        // 末尾から "Z<数字>" を探す
        zPos = -1;
        for (p = lengthOf(basename) - 1; p >= 0; p--) {
            c = substring(basename, p, p + 1);
            if (c == "Z") {
                zPos = p;
                break;
            }
        }
        if (zPos == -1) {
            print("スキップ (Zが見つかりません): " + fileName);
            nSkipped++;
            continue;
        }
        zStr = substring(basename, zPos + 1);
        if (!matches(zStr, "[0-9]+")) {
            print("スキップ (Z番号が数値ではありません): " + fileName);
            nSkipped++;
            continue;
        }
        zIndex = parseInt(zStr);

        // "Z<数字>" の直前の "_" を探して、そこまでを元ファイル名とする
        cPos = -1;
        for (p = zPos - 1; p >= 0; p--) {
            c = substring(basename, p, p + 1);
            if (c == "_") {
                cPos = p;
                break;
            }
        }
        if (cPos == -1) {
            print("スキップ (区切り文字が見つかりません): " + fileName);
            nSkipped++;
            continue;
        }
        oibBase = substring(basename, 0, cPos);

        // tif作成時にだけ付けた波長サフィックスを、元ファイル名から取り除く
        // (表記ゆれがあり得るので、候補を順番に試す)
        for (s = 0; s < wavelengthSuffixCandidates.length; s++) {
            suffix = wavelengthSuffixCandidates[s];
            if (suffix != "" && endsWith(oibBase, suffix)) {
                oibBase = substring(oibBase, 0, lengthOf(oibBase) - lengthOf(suffix));
                s = wavelengthSuffixCandidates.length; // 一致したらループを抜ける
            }
        }

        // olympusDir直下だけでなく、サブフォルダの中も再帰的に探す
        // (.oif は本体が小さいメタデータファイルで、実データは同名の
        //  "○○○.oif.files" フォルダ内に入っているが、Bio-Formatsは
        //  .oifファイルを指定するだけで中身を自動で読んでくれる)
        oibExtensions = newArray(".oib", ".oir", ".oif");
        oibPath = "";
        for (e = 0; e < oibExtensions.length; e++) {
            oibPath = findFileRecursive(olympusDir, oibBase + oibExtensions[e]);
            if (oibPath != "") break;
        }
        if (oibPath == "") {
            print("スキップ (元のOlympusファイルが見つかりません): " + oibBase);
            nSkipped++;
            continue;
        }

        outPath = outputDir + basename + ".tif";

        run("Bio-Formats Importer",
            "open=[" + oibPath + "] color_mode=Default " +
            "specify_range " +
            "c_begin=" + targetChannel + " c_end=" + targetChannel + " c_step=1 " +
            "z_begin=" + zIndex + " z_end=" + zIndex + " z_step=1 " +
            "t_begin=1 t_end=1 t_step=1");

        saveAs("Tiff", outPath);
        close();

        print(basename + " -> 16bit書き出し完了 (channel=" + targetChannel + ", z=" + zIndex + ")");
        nDone++;
    }

    showMessage("完了: " + nDone + "件書き出し / " + nSkipped + "件スキップ\n詳細はLogウィンドウを確認してください。");
}

// dir以下をサブフォルダも含めて再帰的に探索し、
// ファイル名がtargetNameと一致する最初のファイルのフルパスを返す。
// 見つからない場合は空文字を返す。
function findFileRecursive(dir, targetName) {
    list = getFileList(dir);
    for (idx = 0; idx < list.length; idx++) {
        entry = list[idx];
        path = dir + entry;
        if (endsWith(entry, "/")) {
            found = findFileRecursive(path, targetName);
            if (found != "") return found;
        } else if (entry == targetName) {
            return path;
        }
    }
    return "";
}

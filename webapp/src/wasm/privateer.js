var privateer_module = (() => {
  var _scriptDir = import.meta.url;

  return async function (moduleArg = {}) {
    var Module = moduleArg;
    var readyPromiseResolve, readyPromiseReject;
    Module["ready"] = new Promise((resolve, reject) => {
      readyPromiseResolve = resolve;
      readyPromiseReject = reject;
    });
    [
      "_main",
      "getExceptionMessage",
      "___get_exception_message",
      "_free",
      "___cpp_exception",
      "___cxa_increment_exception_refcount",
      "___cxa_decrement_exception_refcount",
      "___thrown_object_from_unwind_exception",
      "_fflush",
      "__embind_initialize_bindings",
      "___set_stack_limits",
      "onRuntimeInitialized",
    ].forEach((prop) => {
      if (!Object.getOwnPropertyDescriptor(Module["ready"], prop)) {
        Object.defineProperty(Module["ready"], prop, {
          get: () =>
            abort(
              "You are getting " +
                prop +
                " on the Promise object, instead of the instance. Use .then() to get called back with the instance, see the MODULARIZE docs in src/settings.js",
            ),
          set: () =>
            abort(
              "You are setting " +
                prop +
                " on the Promise object, instead of the instance. Use .then() to get called back with the instance, see the MODULARIZE docs in src/settings.js",
            ),
        });
      }
    });
    if (!Module.expectedDataFileDownloads) {
      Module.expectedDataFileDownloads = 0;
    }
    Module.expectedDataFileDownloads++;
    (function () {
      if (Module["ENVIRONMENT_IS_PTHREAD"] || Module["$ww"]) return;
      var loadPackage = function (metadata) {
        var PACKAGE_PATH = "";
        if (typeof window === "object") {
          PACKAGE_PATH = window["encodeURIComponent"](
            window.location.pathname
              .toString()
              .substring(
                0,
                window.location.pathname.toString().lastIndexOf("/"),
              ) + "/",
          );
        } else if (
          typeof process === "undefined" &&
          typeof location !== "undefined"
        ) {
          PACKAGE_PATH = encodeURIComponent(
            location.pathname
              .toString()
              .substring(0, location.pathname.toString().lastIndexOf("/")) +
              "/",
          );
        }
        var PACKAGE_NAME = "privateer.data";
        var REMOTE_PACKAGE_BASE = "privateer.data";
        if (
          typeof Module["locateFilePackage"] === "function" &&
          !Module["locateFile"]
        ) {
          Module["locateFile"] = Module["locateFilePackage"];
          err(
            "warning: you defined Module.locateFilePackage, that has been renamed to Module.locateFile (using your locateFilePackage for now)",
          );
        }
        var REMOTE_PACKAGE_NAME = Module["locateFile"]
          ? Module["locateFile"](REMOTE_PACKAGE_BASE, "")
          : REMOTE_PACKAGE_BASE;
        var REMOTE_PACKAGE_SIZE = metadata["remote_package_size"];
        function fetchRemotePackage(
          packageName,
          packageSize,
          callback,
          errback,
        ) {
          if (
            typeof process === "object" &&
            typeof process.versions === "object" &&
            typeof process.versions.node === "string"
          ) {
            require("fs").readFile(packageName, function (err, contents) {
              if (err) {
                errback(err);
              } else {
                callback(contents.buffer);
              }
            });
            return;
          }
          var xhr = new XMLHttpRequest();
          xhr.open("GET", packageName, true);
          xhr.responseType = "arraybuffer";
          xhr.onprogress = function (event) {
            var url = packageName;
            var size = packageSize;
            if (event.total) size = event.total;
            if (event.loaded) {
              if (!xhr.addedTotal) {
                xhr.addedTotal = true;
                if (!Module.dataFileDownloads) Module.dataFileDownloads = {};
                Module.dataFileDownloads[url] = {
                  loaded: event.loaded,
                  total: size,
                };
              } else {
                Module.dataFileDownloads[url].loaded = event.loaded;
              }
              var total = 0;
              var loaded = 0;
              var num = 0;
              for (var download in Module.dataFileDownloads) {
                var data = Module.dataFileDownloads[download];
                total += data.total;
                loaded += data.loaded;
                num++;
              }
              total = Math.ceil(
                (total * Module.expectedDataFileDownloads) / num,
              );
              if (Module["setStatus"])
                Module["setStatus"](`Downloading data... (${loaded}/${total})`);
            } else if (!Module.dataFileDownloads) {
              if (Module["setStatus"])
                Module["setStatus"]("Downloading data...");
            }
          };
          xhr.onerror = function (event) {
            throw new Error("NetworkError for: " + packageName);
          };
          xhr.onload = function (event) {
            if (
              xhr.status == 200 ||
              xhr.status == 304 ||
              xhr.status == 206 ||
              (xhr.status == 0 && xhr.response)
            ) {
              var packageData = xhr.response;
              callback(packageData);
            } else {
              throw new Error(xhr.statusText + " : " + xhr.responseURL);
            }
          };
          xhr.send(null);
        }
        function handleError(error) {
          console.error("package error:", error);
        }
        var fetchedCallback = null;
        var fetched = Module["getPreloadedPackage"]
          ? Module["getPreloadedPackage"](
              REMOTE_PACKAGE_NAME,
              REMOTE_PACKAGE_SIZE,
            )
          : null;
        if (!fetched)
          fetchRemotePackage(
            REMOTE_PACKAGE_NAME,
            REMOTE_PACKAGE_SIZE,
            function (data) {
              if (fetchedCallback) {
                fetchedCallback(data);
                fetchedCallback = null;
              } else {
                fetched = data;
              }
            },
            handleError,
          );
        function runWithFS() {
          function assert(check, msg) {
            if (!check) throw msg + new Error().stack;
          }
          Module["FS_createPath"]("/", "unprocessed_files", true, true);
          function DataRequest(start, end, audio) {
            this.start = start;
            this.end = end;
            this.audio = audio;
          }
          DataRequest.prototype = {
            requests: {},
            open: function (mode, name) {
              this.name = name;
              this.requests[name] = this;
              Module["addRunDependency"](`fp ${this.name}`);
            },
            send: function () {},
            onload: function () {
              var byteArray = this.byteArray.subarray(this.start, this.end);
              this.finish(byteArray);
            },
            finish: function (byteArray) {
              var that = this;
              Module["FS_createDataFile"](
                this.name,
                null,
                byteArray,
                true,
                true,
                true,
              );
              Module["removeRunDependency"](`fp ${that.name}`);
              this.requests[this.name] = null;
            },
          };
          var files = metadata["files"];
          for (var i = 0; i < files.length; ++i) {
            new DataRequest(
              files[i]["start"],
              files[i]["end"],
              files[i]["audio"] || 0,
            ).open("GET", files[i]["filename"]);
          }
          function processPackageData(arrayBuffer) {
            assert(arrayBuffer, "Loading data file failed.");
            assert(
              arrayBuffer.constructor.name === ArrayBuffer.name,
              "bad input to processPackageData",
            );
            var byteArray = new Uint8Array(arrayBuffer);
            DataRequest.prototype.byteArray = byteArray;
            var files = metadata["files"];
            for (var i = 0; i < files.length; ++i) {
              DataRequest.prototype.requests[files[i].filename].onload();
            }
            Module["removeRunDependency"]("datafile_privateer.data");
          }
          Module["addRunDependency"]("datafile_privateer.data");
          if (!Module.preloadResults) Module.preloadResults = {};
          Module.preloadResults[PACKAGE_NAME] = { fromCache: false };
          if (fetched) {
            processPackageData(fetched);
            fetched = null;
          } else {
            fetchedCallback = processPackageData;
          }
        }
        if (Module["calledRun"]) {
          runWithFS();
        } else {
          if (!Module["preRun"]) Module["preRun"] = [];
          Module["preRun"].push(runWithFS);
        }
      };
      loadPackage({
        files: [
          {
            filename: "/privateer_torsion_database.json",
            start: 0,
            end: 2313860,
          },
          {
            filename: "/privateer_torsion_statistics.json",
            start: 2313860,
            end: 2326506,
          },
          {
            filename: "/privateer_torsions_z_score_database.json",
            start: 2326506,
            end: 10108540,
          },
          {
            filename: "/unprocessed_files/ASN-NAG_reduced.json",
            start: 10108540,
            end: 10757358,
          },
          {
            filename: "/unprocessed_files/ASN-NGA_reduced.json",
            start: 10757358,
            end: 10758631,
          },
          {
            filename: "/unprocessed_files/BGC-BGC_reduced.json",
            start: 10758631,
            end: 10760275,
          },
          {
            filename: "/unprocessed_files/BGC-GAL_reduced.json",
            start: 10760275,
            end: 10760795,
          },
          {
            filename: "/unprocessed_files/BGC-GLA_reduced.json",
            start: 10760795,
            end: 10761054,
          },
          {
            filename: "/unprocessed_files/BGC-XYS_reduced.json",
            start: 10761054,
            end: 10761397,
          },
          {
            filename: "/unprocessed_files/BMA-BMA_reduced.json",
            start: 10761397,
            end: 10762693,
          },
          {
            filename: "/unprocessed_files/BMA-GAL_reduced.json",
            start: 10762693,
            end: 10762868,
          },
          {
            filename: "/unprocessed_files/BMA-MAN_reduced.json",
            start: 10762868,
            end: 10888263,
          },
          {
            filename: "/unprocessed_files/BMA-NAG_reduced.json",
            start: 10888263,
            end: 10888862,
          },
          {
            filename: "/unprocessed_files/BMA-NDG_reduced.json",
            start: 10888862,
            end: 10888948,
          },
          {
            filename: "/unprocessed_files/BMA-XYP_reduced.json",
            start: 10888948,
            end: 10889890,
          },
          {
            filename: "/unprocessed_files/FRU-FRU_reduced.json",
            start: 10889890,
            end: 10890666,
          },
          {
            filename: "/unprocessed_files/FRU-GLC_reduced.json",
            start: 10890666,
            end: 10891785,
          },
          {
            filename: "/unprocessed_files/FUC-BGC_reduced.json",
            start: 10891785,
            end: 10892733,
          },
          {
            filename: "/unprocessed_files/GAL-A2G_reduced.json",
            start: 10892733,
            end: 10892991,
          },
          {
            filename: "/unprocessed_files/GAL-FUC_reduced.json",
            start: 10892991,
            end: 10893162,
          },
          {
            filename: "/unprocessed_files/GAL-GAL_reduced.json",
            start: 10893162,
            end: 10894023,
          },
          {
            filename: "/unprocessed_files/GAL-GLA_reduced.json",
            start: 10894023,
            end: 10894532,
          },
          {
            filename: "/unprocessed_files/GAL-NAG_reduced.json",
            start: 10894532,
            end: 10894961,
          },
          {
            filename: "/unprocessed_files/GAL-SIA_reduced.json",
            start: 10894961,
            end: 10899540,
          },
          {
            filename: "/unprocessed_files/GLA-A2G_reduced.json",
            start: 10899540,
            end: 10899714,
          },
          {
            filename: "/unprocessed_files/GLA-GLA_reduced.json",
            start: 10899714,
            end: 10900142,
          },
          {
            filename: "/unprocessed_files/GLC-BGC_reduced.json",
            start: 10900142,
            end: 10900319,
          },
          {
            filename: "/unprocessed_files/GLC-FRU_reduced.json",
            start: 10900319,
            end: 10900493,
          },
          {
            filename: "/unprocessed_files/GLC-GAL_reduced.json",
            start: 10900493,
            end: 10900581,
          },
          {
            filename: "/unprocessed_files/GLC-GLA_reduced.json",
            start: 10900581,
            end: 10900752,
          },
          {
            filename: "/unprocessed_files/GLC-GLC_reduced.json",
            start: 10900752,
            end: 10901525,
          },
          {
            filename: "/unprocessed_files/GLC-MAN_reduced.json",
            start: 10901525,
            end: 10901612,
          },
          {
            filename: "/unprocessed_files/GLC-XYS_reduced.json",
            start: 10901612,
            end: 10901699,
          },
          {
            filename: "/unprocessed_files/MAG-GAL_reduced.json",
            start: 10901699,
            end: 10901787,
          },
          {
            filename: "/unprocessed_files/MAN-BMA_reduced.json",
            start: 10901787,
            end: 10901875,
          },
          {
            filename: "/unprocessed_files/MAN-GLC_reduced.json",
            start: 10901875,
            end: 10902549,
          },
          {
            filename: "/unprocessed_files/MAN-MAN_reduced.json",
            start: 10902549,
            end: 10987921,
          },
          {
            filename: "/unprocessed_files/MAN-NAG_reduced.json",
            start: 10987921,
            end: 11008593,
          },
          {
            filename: "/unprocessed_files/NAG-BMA_reduced.json",
            start: 11008593,
            end: 11152243,
          },
          {
            filename: "/unprocessed_files/NAG-FCA_reduced.json",
            start: 11152243,
            end: 11152331,
          },
          {
            filename: "/unprocessed_files/NAG-FUC_reduced.json",
            start: 11152331,
            end: 11206617,
          },
          {
            filename: "/unprocessed_files/NAG-FUL_reduced.json",
            start: 11206617,
            end: 11208669,
          },
          {
            filename: "/unprocessed_files/NAG-GAL_reduced.json",
            start: 11208669,
            end: 11215929,
          },
          {
            filename: "/unprocessed_files/NAG-GLC_reduced.json",
            start: 11215929,
            end: 11216018,
          },
          {
            filename: "/unprocessed_files/NAG-MAN_reduced.json",
            start: 11216018,
            end: 11216364,
          },
          {
            filename: "/unprocessed_files/NAG-NAG_reduced.json",
            start: 11216364,
            end: 11545985,
          },
          {
            filename: "/unprocessed_files/NAG-NDG_reduced.json",
            start: 11545985,
            end: 11546494,
          },
          {
            filename: "/unprocessed_files/NAG-SIA_reduced.json",
            start: 11546494,
            end: 11546667,
          },
          {
            filename: "/unprocessed_files/NDG-BMA_reduced.json",
            start: 11546667,
            end: 11546928,
          },
          {
            filename: "/unprocessed_files/NDG-GAL_reduced.json",
            start: 11546928,
            end: 11547102,
          },
          {
            filename: "/unprocessed_files/NDG-NAG_reduced.json",
            start: 11547102,
            end: 11547278,
          },
          {
            filename: "/unprocessed_files/NGA-GAL_reduced.json",
            start: 11547278,
            end: 11547454,
          },
          {
            filename: "/unprocessed_files/NGA-NAG_reduced.json",
            start: 11547454,
            end: 11547542,
          },
          {
            filename: "/unprocessed_files/PA1-GCS_reduced.json",
            start: 11547542,
            end: 11547631,
          },
          {
            filename: "/unprocessed_files/XYP-GCV_reduced.json",
            start: 11547631,
            end: 11547719,
          },
          {
            filename: "/unprocessed_files/XYP-XYP_reduced.json",
            start: 11547719,
            end: 11548065,
          },
        ],
        remote_package_size: 11548065,
      });
    })();
    if (Module["ENVIRONMENT_IS_PTHREAD"] || Module["$ww"])
      Module["preRun"] = [];
    var necessaryPreJSTasks = Module["preRun"].slice();
    if (!Module["preRun"])
      throw "Module.preRun should exist because file support used it; did a pre-js delete it?";
    necessaryPreJSTasks.forEach(function (task) {
      if (Module["preRun"].indexOf(task) < 0)
        throw "All preRun tasks that exist before user pre-js code should remain after; did you replace Module or modify Module.preRun?";
    });
    var moduleOverrides = Object.assign({}, Module);
    var arguments_ = [];
    var thisProgram = "./this.program";
    var quit_ = (status, toThrow) => {
      throw toThrow;
    };
    var ENVIRONMENT_IS_WEB = typeof window == "object";
    var ENVIRONMENT_IS_WORKER = typeof importScripts == "function";
    var ENVIRONMENT_IS_NODE =
      typeof process == "object" &&
      typeof process.versions == "object" &&
      typeof process.versions.node == "string";
    var ENVIRONMENT_IS_SHELL =
      !ENVIRONMENT_IS_WEB && !ENVIRONMENT_IS_NODE && !ENVIRONMENT_IS_WORKER;
    if (Module["ENVIRONMENT"]) {
      throw new Error(
        "Module.ENVIRONMENT has been deprecated. To force the environment, use the ENVIRONMENT compile-time option (for example, -sENVIRONMENT=web or -sENVIRONMENT=node)",
      );
    }
    var scriptDirectory = "";
    function locateFile(path) {
      if (Module["locateFile"]) {
        return Module["locateFile"](path, scriptDirectory);
      }
      return scriptDirectory + path;
    }
    var read_, readAsync, readBinary, setWindowTitle;
    if (ENVIRONMENT_IS_NODE) {
      if (
        typeof process == "undefined" ||
        !process.release ||
        process.release.name !== "node"
      )
        throw new Error(
          "not compiled for this environment (did you build to HTML and try to run it not on the web, or set ENVIRONMENT to something - like node - and run it someplace else - like on the web?)",
        );
      var nodeVersion = process.versions.node;
      var numericVersion = nodeVersion.split(".").slice(0, 3);
      numericVersion =
        numericVersion[0] * 1e4 +
        numericVersion[1] * 100 +
        numericVersion[2].split("-")[0] * 1;
      if (numericVersion < 16e4) {
        throw new Error(
          "This emscripten-generated code requires node v16.0.0 (detected v" +
            nodeVersion +
            ")",
        );
      }
      const { createRequire: createRequire } = await import("module");
      var require = createRequire(import.meta.url);
      var fs = require("fs");
      var nodePath = require("path");
      if (ENVIRONMENT_IS_WORKER) {
        scriptDirectory = nodePath.dirname(scriptDirectory) + "/";
      } else {
        scriptDirectory = require("url").fileURLToPath(
          new URL("./", import.meta.url),
        );
      }
      read_ = (filename, binary) => {
        filename = isFileURI(filename)
          ? new URL(filename)
          : nodePath.normalize(filename);
        return fs.readFileSync(filename, binary ? undefined : "utf8");
      };
      readBinary = (filename) => {
        var ret = read_(filename, true);
        if (!ret.buffer) {
          ret = new Uint8Array(ret);
        }
        assert(ret.buffer);
        return ret;
      };
      readAsync = (filename, onload, onerror, binary = true) => {
        filename = isFileURI(filename)
          ? new URL(filename)
          : nodePath.normalize(filename);
        fs.readFile(filename, binary ? undefined : "utf8", (err, data) => {
          if (err) onerror(err);
          else onload(binary ? data.buffer : data);
        });
      };
      if (!Module["thisProgram"] && process.argv.length > 1) {
        thisProgram = process.argv[1].replace(/\\/g, "/");
      }
      arguments_ = process.argv.slice(2);
      quit_ = (status, toThrow) => {
        process.exitCode = status;
        throw toThrow;
      };
      Module["inspect"] = () => "[Emscripten Module object]";
    } else if (ENVIRONMENT_IS_SHELL) {
      if (
        (typeof process == "object" && typeof require === "function") ||
        typeof window == "object" ||
        typeof importScripts == "function"
      )
        throw new Error(
          "not compiled for this environment (did you build to HTML and try to run it not on the web, or set ENVIRONMENT to something - like node - and run it someplace else - like on the web?)",
        );
      if (typeof read != "undefined") {
        read_ = (f) => read(f);
      }
      readBinary = (f) => {
        let data;
        if (typeof readbuffer == "function") {
          return new Uint8Array(readbuffer(f));
        }
        data = read(f, "binary");
        assert(typeof data == "object");
        return data;
      };
      readAsync = (f, onload, onerror) => {
        setTimeout(() => onload(readBinary(f)));
      };
      if (typeof clearTimeout == "undefined") {
        globalThis.clearTimeout = (id) => {};
      }
      if (typeof setTimeout == "undefined") {
        globalThis.setTimeout = (f) => (typeof f == "function" ? f() : abort());
      }
      if (typeof scriptArgs != "undefined") {
        arguments_ = scriptArgs;
      } else if (typeof arguments != "undefined") {
        arguments_ = arguments;
      }
      if (typeof quit == "function") {
        quit_ = (status, toThrow) => {
          setTimeout(() => {
            if (!(toThrow instanceof ExitStatus)) {
              let toLog = toThrow;
              if (toThrow && typeof toThrow == "object" && toThrow.stack) {
                toLog = [toThrow, toThrow.stack];
              }
              err(`exiting due to exception: ${toLog}`);
            }
            quit(status);
          });
          throw toThrow;
        };
      }
      if (typeof print != "undefined") {
        if (typeof console == "undefined") console = {};
        console.log = print;
        console.warn = console.error =
          typeof printErr != "undefined" ? printErr : print;
      }
    } else if (ENVIRONMENT_IS_WEB || ENVIRONMENT_IS_WORKER) {
      if (ENVIRONMENT_IS_WORKER) {
        scriptDirectory = self.location.href;
      } else if (typeof document != "undefined" && document.currentScript) {
        scriptDirectory = document.currentScript.src;
      }
      if (_scriptDir) {
        scriptDirectory = _scriptDir;
      }
      if (scriptDirectory.indexOf("blob:") !== 0) {
        scriptDirectory = scriptDirectory.substr(
          0,
          scriptDirectory.replace(/[?#].*/, "").lastIndexOf("/") + 1,
        );
      } else {
        scriptDirectory = "";
      }
      if (!(typeof window == "object" || typeof importScripts == "function"))
        throw new Error(
          "not compiled for this environment (did you build to HTML and try to run it not on the web, or set ENVIRONMENT to something - like node - and run it someplace else - like on the web?)",
        );
      {
        read_ = (url) => {
          var xhr = new XMLHttpRequest();
          xhr.open("GET", url, false);
          xhr.send(null);
          return xhr.responseText;
        };
        if (ENVIRONMENT_IS_WORKER) {
          readBinary = (url) => {
            var xhr = new XMLHttpRequest();
            xhr.open("GET", url, false);
            xhr.responseType = "arraybuffer";
            xhr.send(null);
            return new Uint8Array(xhr.response);
          };
        }
        readAsync = (url, onload, onerror) => {
          var xhr = new XMLHttpRequest();
          xhr.open("GET", url, true);
          xhr.responseType = "arraybuffer";
          xhr.onload = () => {
            if (xhr.status == 200 || (xhr.status == 0 && xhr.response)) {
              onload(xhr.response);
              return;
            }
            onerror();
          };
          xhr.onerror = onerror;
          xhr.send(null);
        };
      }
      setWindowTitle = (title) => (document.title = title);
    } else {
      throw new Error("environment detection error");
    }
    var out = Module["print"] || console.log.bind(console);
    var err = Module["printErr"] || console.error.bind(console);
    Object.assign(Module, moduleOverrides);
    moduleOverrides = null;
    checkIncomingModuleAPI();
    if (Module["arguments"]) arguments_ = Module["arguments"];
    legacyModuleProp("arguments", "arguments_");
    if (Module["thisProgram"]) thisProgram = Module["thisProgram"];
    legacyModuleProp("thisProgram", "thisProgram");
    if (Module["quit"]) quit_ = Module["quit"];
    legacyModuleProp("quit", "quit_");
    assert(
      typeof Module["memoryInitializerPrefixURL"] == "undefined",
      "Module.memoryInitializerPrefixURL option was removed, use Module.locateFile instead",
    );
    assert(
      typeof Module["pthreadMainPrefixURL"] == "undefined",
      "Module.pthreadMainPrefixURL option was removed, use Module.locateFile instead",
    );
    assert(
      typeof Module["cdInitializerPrefixURL"] == "undefined",
      "Module.cdInitializerPrefixURL option was removed, use Module.locateFile instead",
    );
    assert(
      typeof Module["filePackagePrefixURL"] == "undefined",
      "Module.filePackagePrefixURL option was removed, use Module.locateFile instead",
    );
    assert(
      typeof Module["read"] == "undefined",
      "Module.read option was removed (modify read_ in JS)",
    );
    assert(
      typeof Module["readAsync"] == "undefined",
      "Module.readAsync option was removed (modify readAsync in JS)",
    );
    assert(
      typeof Module["readBinary"] == "undefined",
      "Module.readBinary option was removed (modify readBinary in JS)",
    );
    assert(
      typeof Module["setWindowTitle"] == "undefined",
      "Module.setWindowTitle option was removed (modify setWindowTitle in JS)",
    );
    assert(
      typeof Module["TOTAL_MEMORY"] == "undefined",
      "Module.TOTAL_MEMORY has been renamed Module.INITIAL_MEMORY",
    );
    legacyModuleProp("read", "read_");
    legacyModuleProp("readAsync", "readAsync");
    legacyModuleProp("readBinary", "readBinary");
    legacyModuleProp("setWindowTitle", "setWindowTitle");
    assert(
      !ENVIRONMENT_IS_SHELL,
      "shell environment detected but not enabled at build time.  Add 'shell' to `-sENVIRONMENT` to enable.",
    );
    var wasmBinary;
    if (Module["wasmBinary"]) wasmBinary = Module["wasmBinary"];
    legacyModuleProp("wasmBinary", "wasmBinary");
    var noExitRuntime = Module["noExitRuntime"] || true;
    legacyModuleProp("noExitRuntime", "noExitRuntime");
    if (typeof WebAssembly != "object") {
      abort("no native wasm support detected");
    }
    var wasmMemory;
    var ABORT = false;
    var EXITSTATUS;
    function assert(condition, text) {
      if (!condition) {
        abort("Assertion failed" + (text ? ": " + text : ""));
      }
    }
    var HEAP8, HEAPU8, HEAP16, HEAPU16, HEAP32, HEAPU32, HEAPF32, HEAPF64;
    function updateMemoryViews() {
      var b = wasmMemory.buffer;
      Module["HEAP8"] = HEAP8 = new Int8Array(b);
      Module["HEAP16"] = HEAP16 = new Int16Array(b);
      Module["HEAP32"] = HEAP32 = new Int32Array(b);
      Module["HEAPU8"] = HEAPU8 = new Uint8Array(b);
      Module["HEAPU16"] = HEAPU16 = new Uint16Array(b);
      Module["HEAPU32"] = HEAPU32 = new Uint32Array(b);
      Module["HEAPF32"] = HEAPF32 = new Float32Array(b);
      Module["HEAPF64"] = HEAPF64 = new Float64Array(b);
    }
    assert(
      !Module["STACK_SIZE"],
      "STACK_SIZE can no longer be set at runtime.  Use -sSTACK_SIZE at link time",
    );
    assert(
      typeof Int32Array != "undefined" &&
        typeof Float64Array !== "undefined" &&
        Int32Array.prototype.subarray != undefined &&
        Int32Array.prototype.set != undefined,
      "JS engine does not provide full typed array support",
    );
    assert(
      !Module["wasmMemory"],
      "Use of `wasmMemory` detected.  Use -sIMPORTED_MEMORY to define wasmMemory externally",
    );
    assert(
      !Module["INITIAL_MEMORY"],
      "Detected runtime INITIAL_MEMORY setting.  Use -sIMPORTED_MEMORY to define wasmMemory dynamically",
    );
    var wasmTable;
    function writeStackCookie() {
      var max = _emscripten_stack_get_end();
      assert((max & 3) == 0);
      if (max == 0) {
        max += 4;
      }
      HEAPU32[max >> 2] = 34821223;
      checkInt32(34821223);
      HEAPU32[(max + 4) >> 2] = 2310721022;
      checkInt32(2310721022);
      HEAPU32[0 >> 2] = 1668509029;
      checkInt32(1668509029);
    }
    function checkStackCookie() {
      if (ABORT) return;
      var max = _emscripten_stack_get_end();
      if (max == 0) {
        max += 4;
      }
      var cookie1 = HEAPU32[max >> 2];
      var cookie2 = HEAPU32[(max + 4) >> 2];
      if (cookie1 != 34821223 || cookie2 != 2310721022) {
        abort(
          `Stack overflow! Stack cookie has been overwritten at ${ptrToString(
            max,
          )}, expected hex dwords 0x89BACDFE and 0x2135467, but received ${ptrToString(
            cookie2,
          )} ${ptrToString(cookie1)}`,
        );
      }
      if (HEAPU32[0 >> 2] != 1668509029) {
        abort(
          "Runtime error: The application has corrupted its heap memory area (address zero)!",
        );
      }
    }
    (function () {
      var h16 = new Int16Array(1);
      var h8 = new Int8Array(h16.buffer);
      h16[0] = 25459;
      if (h8[0] !== 115 || h8[1] !== 99)
        throw "Runtime error: expected the system to be little-endian! (Run with -sSUPPORT_BIG_ENDIAN to bypass)";
    })();
    var __ATPRERUN__ = [];
    var __ATINIT__ = [];
    var __ATPOSTRUN__ = [];
    var runtimeInitialized = false;
    var runtimeKeepaliveCounter = 0;
    function preRun() {
      if (Module["preRun"]) {
        if (typeof Module["preRun"] == "function")
          Module["preRun"] = [Module["preRun"]];
        while (Module["preRun"].length) {
          addOnPreRun(Module["preRun"].shift());
        }
      }
      callRuntimeCallbacks(__ATPRERUN__);
    }
    function initRuntime() {
      assert(!runtimeInitialized);
      runtimeInitialized = true;
      checkStackCookie();
      setStackLimits();
      if (!Module["noFSInit"] && !FS.init.initialized) FS.init();
      FS.ignorePermissions = false;
      TTY.init();
      callRuntimeCallbacks(__ATINIT__);
    }
    function postRun() {
      checkStackCookie();
      if (Module["postRun"]) {
        if (typeof Module["postRun"] == "function")
          Module["postRun"] = [Module["postRun"]];
        while (Module["postRun"].length) {
          addOnPostRun(Module["postRun"].shift());
        }
      }
      callRuntimeCallbacks(__ATPOSTRUN__);
    }
    function addOnPreRun(cb) {
      __ATPRERUN__.unshift(cb);
    }
    function addOnInit(cb) {
      __ATINIT__.unshift(cb);
    }
    function addOnPostRun(cb) {
      __ATPOSTRUN__.unshift(cb);
    }
    assert(
      Math.imul,
      "This browser does not support Math.imul(), build with LEGACY_VM_SUPPORT or POLYFILL_OLD_MATH_FUNCTIONS to add in a polyfill",
    );
    assert(
      Math.fround,
      "This browser does not support Math.fround(), build with LEGACY_VM_SUPPORT or POLYFILL_OLD_MATH_FUNCTIONS to add in a polyfill",
    );
    assert(
      Math.clz32,
      "This browser does not support Math.clz32(), build with LEGACY_VM_SUPPORT or POLYFILL_OLD_MATH_FUNCTIONS to add in a polyfill",
    );
    assert(
      Math.trunc,
      "This browser does not support Math.trunc(), build with LEGACY_VM_SUPPORT or POLYFILL_OLD_MATH_FUNCTIONS to add in a polyfill",
    );
    var runDependencies = 0;
    var runDependencyWatcher = null;
    var dependenciesFulfilled = null;
    var runDependencyTracking = {};
    function getUniqueRunDependency(id) {
      var orig = id;
      while (1) {
        if (!runDependencyTracking[id]) return id;
        id = orig + Math.random();
      }
    }
    function addRunDependency(id) {
      runDependencies++;
      if (Module["monitorRunDependencies"]) {
        Module["monitorRunDependencies"](runDependencies);
      }
      if (id) {
        assert(!runDependencyTracking[id]);
        runDependencyTracking[id] = 1;
        if (
          runDependencyWatcher === null &&
          typeof setInterval != "undefined"
        ) {
          runDependencyWatcher = setInterval(() => {
            if (ABORT) {
              clearInterval(runDependencyWatcher);
              runDependencyWatcher = null;
              return;
            }
            var shown = false;
            for (var dep in runDependencyTracking) {
              if (!shown) {
                shown = true;
                err("still waiting on run dependencies:");
              }
              err("dependency: " + dep);
            }
            if (shown) {
              err("(end of list)");
            }
          }, 1e4);
        }
      } else {
        err("warning: run dependency added without ID");
      }
    }
    function removeRunDependency(id) {
      runDependencies--;
      if (Module["monitorRunDependencies"]) {
        Module["monitorRunDependencies"](runDependencies);
      }
      if (id) {
        assert(runDependencyTracking[id]);
        delete runDependencyTracking[id];
      } else {
        err("warning: run dependency removed without ID");
      }
      if (runDependencies == 0) {
        if (runDependencyWatcher !== null) {
          clearInterval(runDependencyWatcher);
          runDependencyWatcher = null;
        }
        if (dependenciesFulfilled) {
          var callback = dependenciesFulfilled;
          dependenciesFulfilled = null;
          callback();
        }
      }
    }
    function abort(what) {
      if (Module["onAbort"]) {
        Module["onAbort"](what);
      }
      what = "Aborted(" + what + ")";
      err(what);
      ABORT = true;
      EXITSTATUS = 1;
      if (runtimeInitialized) {
        ___trap();
      }
      var e = new WebAssembly.RuntimeError(what);
      readyPromiseReject(e);
      throw e;
    }
    var dataURIPrefix = "data:application/octet-stream;base64,";
    function isDataURI(filename) {
      return filename.startsWith(dataURIPrefix);
    }
    function isFileURI(filename) {
      return filename.startsWith("file://");
    }
    function createExportWrapper(name, fixedasm) {
      return function () {
        var displayName = name;
        var asm = fixedasm;
        if (!fixedasm) {
          asm = Module["asm"];
        }
        assert(
          runtimeInitialized,
          "native function `" +
            displayName +
            "` called before runtime initialization",
        );
        if (!asm[name]) {
          assert(
            asm[name],
            "exported native function `" + displayName + "` not found",
          );
        }
        return asm[name].apply(null, arguments);
      };
    }
    var wasmBinaryFile;
    if (Module["locateFile"]) {
      wasmBinaryFile = "privateer.wasm";
      if (!isDataURI(wasmBinaryFile)) {
        wasmBinaryFile = locateFile(wasmBinaryFile);
      }
    } else {
      wasmBinaryFile = new URL("privateer.wasm", import.meta.url).href;
    }
    function getBinarySync(file) {
      if (file == wasmBinaryFile && wasmBinary) {
        return new Uint8Array(wasmBinary);
      }
      if (readBinary) {
        return readBinary(file);
      }
      throw "both async and sync fetching of the wasm failed";
    }
    function getBinaryPromise(binaryFile) {
      if (!wasmBinary && (ENVIRONMENT_IS_WEB || ENVIRONMENT_IS_WORKER)) {
        if (typeof fetch == "function" && !isFileURI(binaryFile)) {
          return fetch(binaryFile, { credentials: "same-origin" })
            .then((response) => {
              if (!response["ok"]) {
                throw "failed to load wasm binary file at '" + binaryFile + "'";
              }
              return response["arrayBuffer"]();
            })
            .catch(() => getBinarySync(binaryFile));
        } else if (readAsync) {
          return new Promise((resolve, reject) => {
            readAsync(
              binaryFile,
              (response) => resolve(new Uint8Array(response)),
              reject,
            );
          });
        }
      }
      return Promise.resolve().then(() => getBinarySync(binaryFile));
    }
    function instantiateArrayBuffer(binaryFile, imports, receiver) {
      return getBinaryPromise(binaryFile)
        .then((binary) => WebAssembly.instantiate(binary, imports))
        .then((instance) => instance)
        .then(receiver, (reason) => {
          err("failed to asynchronously prepare wasm: " + reason);
          if (isFileURI(wasmBinaryFile)) {
            err(
              "warning: Loading from a file URI (" +
                wasmBinaryFile +
                ") is not supported in most browsers. See https://emscripten.org/docs/getting_started/FAQ.html#how-do-i-run-a-local-webserver-for-testing-why-does-my-program-stall-in-downloading-or-preparing",
            );
          }
          abort(reason);
        });
    }
    function instantiateAsync(binary, binaryFile, imports, callback) {
      if (
        !binary &&
        typeof WebAssembly.instantiateStreaming == "function" &&
        !isDataURI(binaryFile) &&
        !isFileURI(binaryFile) &&
        !ENVIRONMENT_IS_NODE &&
        typeof fetch == "function"
      ) {
        return fetch(binaryFile, { credentials: "same-origin" }).then(
          (response) => {
            var result = WebAssembly.instantiateStreaming(response, imports);
            return result.then(callback, function (reason) {
              err("wasm streaming compile failed: " + reason);
              err("falling back to ArrayBuffer instantiation");
              return instantiateArrayBuffer(binaryFile, imports, callback);
            });
          },
        );
      }
      return instantiateArrayBuffer(binaryFile, imports, callback);
    }
    function createWasm() {
      var info = { env: wasmImports, wasi_snapshot_preview1: wasmImports };
      function receiveInstance(instance, module) {
        var exports = instance.exports;
        Module["asm"] = exports;
        wasmMemory = Module["asm"]["memory"];
        assert(wasmMemory, "memory not found in wasm exports");
        updateMemoryViews();
        wasmTable = Module["asm"]["__indirect_function_table"];
        assert(wasmTable, "table not found in wasm exports");
        addOnInit(Module["asm"]["__wasm_call_ctors"]);
        removeRunDependency("wasm-instantiate");
        return exports;
      }
      addRunDependency("wasm-instantiate");
      var trueModule = Module;
      function receiveInstantiationResult(result) {
        assert(
          Module === trueModule,
          "the Module object should not be replaced during async compilation - perhaps the order of HTML elements is wrong?",
        );
        trueModule = null;
        receiveInstance(result["instance"]);
      }
      if (Module["instantiateWasm"]) {
        try {
          return Module["instantiateWasm"](info, receiveInstance);
        } catch (e) {
          err("Module.instantiateWasm callback failed with error: " + e);
          readyPromiseReject(e);
        }
      }
      instantiateAsync(
        wasmBinary,
        wasmBinaryFile,
        info,
        receiveInstantiationResult,
      ).catch(readyPromiseReject);
      return {};
    }
    var tempDouble;
    var tempI64;
    function legacyModuleProp(prop, newName) {
      if (!Object.getOwnPropertyDescriptor(Module, prop)) {
        Object.defineProperty(Module, prop, {
          configurable: true,
          get() {
            abort(
              "Module." +
                prop +
                " has been replaced with plain " +
                newName +
                " (the initial value can be provided on Module, but after startup the value is only looked for on a local variable of that name)",
            );
          },
        });
      }
    }
    function ignoredModuleProp(prop) {
      if (Object.getOwnPropertyDescriptor(Module, prop)) {
        abort(
          "`Module." +
            prop +
            "` was supplied but `" +
            prop +
            "` not included in INCOMING_MODULE_JS_API",
        );
      }
    }
    function isExportedByForceFilesystem(name) {
      return (
        name === "FS_createPath" ||
        name === "FS_createDataFile" ||
        name === "FS_createPreloadedFile" ||
        name === "FS_unlink" ||
        name === "addRunDependency" ||
        name === "FS_createLazyFile" ||
        name === "FS_createDevice" ||
        name === "removeRunDependency"
      );
    }
    function missingGlobal(sym, msg) {
      if (typeof globalThis !== "undefined") {
        Object.defineProperty(globalThis, sym, {
          configurable: true,
          get() {
            warnOnce(
              "`" + sym + "` is not longer defined by emscripten. " + msg,
            );
            return undefined;
          },
        });
      }
    }
    missingGlobal("buffer", "Please use HEAP8.buffer or wasmMemory.buffer");
    function missingLibrarySymbol(sym) {
      if (
        typeof globalThis !== "undefined" &&
        !Object.getOwnPropertyDescriptor(globalThis, sym)
      ) {
        Object.defineProperty(globalThis, sym, {
          configurable: true,
          get() {
            var msg =
              "`" +
              sym +
              "` is a library symbol and not included by default; add it to your library.js __deps or to DEFAULT_LIBRARY_FUNCS_TO_INCLUDE on the command line";
            var librarySymbol = sym;
            if (!librarySymbol.startsWith("_")) {
              librarySymbol = "$" + sym;
            }
            msg +=
              " (e.g. -sDEFAULT_LIBRARY_FUNCS_TO_INCLUDE='" +
              librarySymbol +
              "')";
            if (isExportedByForceFilesystem(sym)) {
              msg +=
                ". Alternatively, forcing filesystem support (-sFORCE_FILESYSTEM) can export this for you";
            }
            warnOnce(msg);
            return undefined;
          },
        });
      }
      unexportedRuntimeSymbol(sym);
    }
    function unexportedRuntimeSymbol(sym) {
      if (!Object.getOwnPropertyDescriptor(Module, sym)) {
        Object.defineProperty(Module, sym, {
          configurable: true,
          get() {
            var msg =
              "'" +
              sym +
              "' was not exported. add it to EXPORTED_RUNTIME_METHODS (see the Emscripten FAQ)";
            if (isExportedByForceFilesystem(sym)) {
              msg +=
                ". Alternatively, forcing filesystem support (-sFORCE_FILESYSTEM) can export this for you";
            }
            abort(msg);
          },
        });
      }
    }
    var MAX_UINT8 = 2 ** 8 - 1;
    var MAX_UINT16 = 2 ** 16 - 1;
    var MAX_UINT32 = 2 ** 32 - 1;
    var MAX_UINT53 = 2 ** 53 - 1;
    var MAX_UINT64 = 2 ** 64 - 1;
    var MIN_INT8 = -(2 ** (8 - 1)) + 1;
    var MIN_INT16 = -(2 ** (16 - 1)) + 1;
    var MIN_INT32 = -(2 ** (32 - 1)) + 1;
    var MIN_INT53 = -(2 ** (53 - 1)) + 1;
    var MIN_INT64 = -(2 ** (64 - 1)) + 1;
    function checkInt(value, bits, min, max) {
      assert(
        Number.isInteger(Number(value)),
        "attempt to write non-integer (" + value + ") into integer heap",
      );
      assert(
        value <= max,
        "value (" + value + ") too large to write as " + bits + "-bit value",
      );
      assert(
        value >= min,
        "value (" + value + ") too small to write as " + bits + "-bit value",
      );
    }
    var checkInt8 = (value) => checkInt(value, 8, MIN_INT8, MAX_UINT8);
    var checkInt16 = (value) => checkInt(value, 16, MIN_INT16, MAX_UINT16);
    var checkInt32 = (value) => checkInt(value, 32, MIN_INT32, MAX_UINT32);
    var checkInt64 = (value) => checkInt(value, 64, MIN_INT64, MAX_UINT64);
    function ExitStatus(status) {
      this.name = "ExitStatus";
      this.message = `Program terminated with exit(${status})`;
      this.status = status;
    }
    var callRuntimeCallbacks = (callbacks) => {
      while (callbacks.length > 0) {
        callbacks.shift()(Module);
      }
    };
    function getCppExceptionTag() {
      return Module["asm"]["__cpp_exception"];
    }
    function getCppExceptionThrownObjectFromWebAssemblyException(ex) {
      var unwind_header = ex.getArg(getCppExceptionTag(), 0);
      return ___thrown_object_from_unwind_exception(unwind_header);
    }
    var withStackSave = (f) => {
      var stack = stackSave();
      var ret = f();
      stackRestore(stack);
      return ret;
    };
    var lengthBytesUTF8 = (str) => {
      var len = 0;
      for (var i = 0; i < str.length; ++i) {
        var c = str.charCodeAt(i);
        if (c <= 127) {
          len++;
        } else if (c <= 2047) {
          len += 2;
        } else if (c >= 55296 && c <= 57343) {
          len += 4;
          ++i;
        } else {
          len += 3;
        }
      }
      return len;
    };
    var stringToUTF8Array = (str, heap, outIdx, maxBytesToWrite) => {
      assert(typeof str === "string");
      if (!(maxBytesToWrite > 0)) return 0;
      var startIdx = outIdx;
      var endIdx = outIdx + maxBytesToWrite - 1;
      for (var i = 0; i < str.length; ++i) {
        var u = str.charCodeAt(i);
        if (u >= 55296 && u <= 57343) {
          var u1 = str.charCodeAt(++i);
          u = (65536 + ((u & 1023) << 10)) | (u1 & 1023);
        }
        if (u <= 127) {
          if (outIdx >= endIdx) break;
          heap[outIdx++] = u;
        } else if (u <= 2047) {
          if (outIdx + 1 >= endIdx) break;
          heap[outIdx++] = 192 | (u >> 6);
          heap[outIdx++] = 128 | (u & 63);
        } else if (u <= 65535) {
          if (outIdx + 2 >= endIdx) break;
          heap[outIdx++] = 224 | (u >> 12);
          heap[outIdx++] = 128 | ((u >> 6) & 63);
          heap[outIdx++] = 128 | (u & 63);
        } else {
          if (outIdx + 3 >= endIdx) break;
          if (u > 1114111)
            warnOnce(
              "Invalid Unicode code point " +
                ptrToString(u) +
                " encountered when serializing a JS string to a UTF-8 string in wasm memory! (Valid unicode code points should be in range 0-0x10FFFF).",
            );
          heap[outIdx++] = 240 | (u >> 18);
          heap[outIdx++] = 128 | ((u >> 12) & 63);
          heap[outIdx++] = 128 | ((u >> 6) & 63);
          heap[outIdx++] = 128 | (u & 63);
        }
      }
      heap[outIdx] = 0;
      return outIdx - startIdx;
    };
    var stringToUTF8 = (str, outPtr, maxBytesToWrite) => {
      assert(
        typeof maxBytesToWrite == "number",
        "stringToUTF8(str, outPtr, maxBytesToWrite) is missing the third parameter that specifies the length of the output buffer!",
      );
      return stringToUTF8Array(str, HEAPU8, outPtr, maxBytesToWrite);
    };
    var stringToUTF8OnStack = (str) => {
      var size = lengthBytesUTF8(str) + 1;
      var ret = stackAlloc(size);
      stringToUTF8(str, ret, size);
      return ret;
    };
    var UTF8Decoder =
      typeof TextDecoder != "undefined" ? new TextDecoder("utf8") : undefined;
    var UTF8ArrayToString = (heapOrArray, idx, maxBytesToRead) => {
      var endIdx = idx + maxBytesToRead;
      var endPtr = idx;
      while (heapOrArray[endPtr] && !(endPtr >= endIdx)) ++endPtr;
      if (endPtr - idx > 16 && heapOrArray.buffer && UTF8Decoder) {
        return UTF8Decoder.decode(heapOrArray.subarray(idx, endPtr));
      }
      var str = "";
      while (idx < endPtr) {
        var u0 = heapOrArray[idx++];
        if (!(u0 & 128)) {
          str += String.fromCharCode(u0);
          continue;
        }
        var u1 = heapOrArray[idx++] & 63;
        if ((u0 & 224) == 192) {
          str += String.fromCharCode(((u0 & 31) << 6) | u1);
          continue;
        }
        var u2 = heapOrArray[idx++] & 63;
        if ((u0 & 240) == 224) {
          u0 = ((u0 & 15) << 12) | (u1 << 6) | u2;
        } else {
          if ((u0 & 248) != 240)
            warnOnce(
              "Invalid UTF-8 leading byte " +
                ptrToString(u0) +
                " encountered when deserializing a UTF-8 string in wasm memory to a JS string!",
            );
          u0 =
            ((u0 & 7) << 18) |
            (u1 << 12) |
            (u2 << 6) |
            (heapOrArray[idx++] & 63);
        }
        if (u0 < 65536) {
          str += String.fromCharCode(u0);
        } else {
          var ch = u0 - 65536;
          str += String.fromCharCode(55296 | (ch >> 10), 56320 | (ch & 1023));
        }
      }
      return str;
    };
    var UTF8ToString = (ptr, maxBytesToRead) => {
      assert(typeof ptr == "number");
      return ptr ? UTF8ArrayToString(HEAPU8, ptr, maxBytesToRead) : "";
    };
    function demangle(func) {
      demangle.recursionGuard = (demangle.recursionGuard | 0) + 1;
      if (demangle.recursionGuard > 1) return func;
      return withStackSave(function () {
        try {
          var s = func;
          if (s.startsWith("__Z")) s = s.substr(1);
          var buf = stringToUTF8OnStack(s);
          var status = stackAlloc(4);
          var ret = ___cxa_demangle(buf, 0, 0, status);
          if (HEAP32[status >> 2] === 0 && ret) {
            return UTF8ToString(ret);
          }
        } catch (e) {
        } finally {
          _free(ret);
          if (demangle.recursionGuard < 2) --demangle.recursionGuard;
        }
        return func;
      });
    }
    var getExceptionMessageCommon = (ptr) =>
      withStackSave(() => {
        var type_addr_addr = stackAlloc(4);
        var message_addr_addr = stackAlloc(4);
        ___get_exception_message(ptr, type_addr_addr, message_addr_addr);
        var type_addr = HEAPU32[type_addr_addr >> 2];
        var message_addr = HEAPU32[message_addr_addr >> 2];
        var type = UTF8ToString(type_addr);
        _free(type_addr);
        var message;
        if (message_addr) {
          message = UTF8ToString(message_addr);
          _free(message_addr);
        }
        return [type, message];
      });
    function getExceptionMessage(ex) {
      var ptr = getCppExceptionThrownObjectFromWebAssemblyException(ex);
      return getExceptionMessageCommon(ptr);
    }
    Module["getExceptionMessage"] = getExceptionMessage;
    var ptrToString = (ptr) => {
      assert(typeof ptr === "number");
      ptr >>>= 0;
      return "0x" + ptr.toString(16).padStart(8, "0");
    };
    var setStackLimits = () => {
      var stackLow = _emscripten_stack_get_base();
      var stackHigh = _emscripten_stack_get_end();
      ___set_stack_limits(stackLow, stackHigh);
    };
    function jsStackTrace() {
      var error = new Error();
      if (!error.stack) {
        try {
          throw new Error();
        } catch (e) {
          error = e;
        }
        if (!error.stack) {
          return "(no stack trace available)";
        }
      }
      return error.stack.toString();
    }
    function demangleAll(text) {
      var regex = /\b_Z[\w\d_]+/g;
      return text.replace(regex, function (x) {
        var y = demangle(x);
        return x === y ? x : y + " [" + x + "]";
      });
    }
    var warnOnce = (text) => {
      if (!warnOnce.shown) warnOnce.shown = {};
      if (!warnOnce.shown[text]) {
        warnOnce.shown[text] = 1;
        if (ENVIRONMENT_IS_NODE) text = "warning: " + text;
        err(text);
      }
    };
    var ___assert_fail = (condition, filename, line, func) => {
      abort(
        `Assertion failed: ${UTF8ToString(condition)}, at: ` +
          [
            filename ? UTF8ToString(filename) : "unknown filename",
            line,
            func ? UTF8ToString(func) : "unknown function",
          ],
      );
    };
    var ___handle_stack_overflow = (requested) => {
      var base = _emscripten_stack_get_base();
      var end = _emscripten_stack_get_end();
      abort(
        `stack overflow (Attempt to set SP to ${ptrToString(requested)}` +
          `, with stack limits [${ptrToString(end)} - ${ptrToString(base)}` +
          "]). If you require more stack space build with -sSTACK_SIZE=<bytes>",
      );
    };
    var setErrNo = (value) => {
      HEAP32[___errno_location() >> 2] = value;
      checkInt32(value);
      return value;
    };
    var PATH = {
      isAbs: (path) => path.charAt(0) === "/",
      splitPath: (filename) => {
        var splitPathRe =
          /^(\/?|)([\s\S]*?)((?:\.{1,2}|[^\/]+?|)(\.[^.\/]*|))(?:[\/]*)$/;
        return splitPathRe.exec(filename).slice(1);
      },
      normalizeArray: (parts, allowAboveRoot) => {
        var up = 0;
        for (var i = parts.length - 1; i >= 0; i--) {
          var last = parts[i];
          if (last === ".") {
            parts.splice(i, 1);
          } else if (last === "..") {
            parts.splice(i, 1);
            up++;
          } else if (up) {
            parts.splice(i, 1);
            up--;
          }
        }
        if (allowAboveRoot) {
          for (; up; up--) {
            parts.unshift("..");
          }
        }
        return parts;
      },
      normalize: (path) => {
        var isAbsolute = PATH.isAbs(path),
          trailingSlash = path.substr(-1) === "/";
        path = PATH.normalizeArray(
          path.split("/").filter((p) => !!p),
          !isAbsolute,
        ).join("/");
        if (!path && !isAbsolute) {
          path = ".";
        }
        if (path && trailingSlash) {
          path += "/";
        }
        return (isAbsolute ? "/" : "") + path;
      },
      dirname: (path) => {
        var result = PATH.splitPath(path),
          root = result[0],
          dir = result[1];
        if (!root && !dir) {
          return ".";
        }
        if (dir) {
          dir = dir.substr(0, dir.length - 1);
        }
        return root + dir;
      },
      basename: (path) => {
        if (path === "/") return "/";
        path = PATH.normalize(path);
        path = path.replace(/\/$/, "");
        var lastSlash = path.lastIndexOf("/");
        if (lastSlash === -1) return path;
        return path.substr(lastSlash + 1);
      },
      join: function () {
        var paths = Array.prototype.slice.call(arguments);
        return PATH.normalize(paths.join("/"));
      },
      join2: (l, r) => PATH.normalize(l + "/" + r),
    };
    var initRandomFill = () => {
      if (
        typeof crypto == "object" &&
        typeof crypto["getRandomValues"] == "function"
      ) {
        return (view) => crypto.getRandomValues(view);
      } else if (ENVIRONMENT_IS_NODE) {
        try {
          var crypto_module = require("crypto");
          var randomFillSync = crypto_module["randomFillSync"];
          if (randomFillSync) {
            return (view) => crypto_module["randomFillSync"](view);
          }
          var randomBytes = crypto_module["randomBytes"];
          return (view) => (view.set(randomBytes(view.byteLength)), view);
        } catch (e) {}
      }
      abort(
        "no cryptographic support found for randomDevice. consider polyfilling it if you want to use something insecure like Math.random(), e.g. put this in a --pre-js: var crypto = { getRandomValues: (array) => { for (var i = 0; i < array.length; i++) array[i] = (Math.random()*256)|0 } };",
      );
    };
    var randomFill = (view) => (randomFill = initRandomFill())(view);
    var PATH_FS = {
      resolve: function () {
        var resolvedPath = "",
          resolvedAbsolute = false;
        for (var i = arguments.length - 1; i >= -1 && !resolvedAbsolute; i--) {
          var path = i >= 0 ? arguments[i] : FS.cwd();
          if (typeof path != "string") {
            throw new TypeError("Arguments to path.resolve must be strings");
          } else if (!path) {
            return "";
          }
          resolvedPath = path + "/" + resolvedPath;
          resolvedAbsolute = PATH.isAbs(path);
        }
        resolvedPath = PATH.normalizeArray(
          resolvedPath.split("/").filter((p) => !!p),
          !resolvedAbsolute,
        ).join("/");
        return (resolvedAbsolute ? "/" : "") + resolvedPath || ".";
      },
      relative: (from, to) => {
        from = PATH_FS.resolve(from).substr(1);
        to = PATH_FS.resolve(to).substr(1);
        function trim(arr) {
          var start = 0;
          for (; start < arr.length; start++) {
            if (arr[start] !== "") break;
          }
          var end = arr.length - 1;
          for (; end >= 0; end--) {
            if (arr[end] !== "") break;
          }
          if (start > end) return [];
          return arr.slice(start, end - start + 1);
        }
        var fromParts = trim(from.split("/"));
        var toParts = trim(to.split("/"));
        var length = Math.min(fromParts.length, toParts.length);
        var samePartsLength = length;
        for (var i = 0; i < length; i++) {
          if (fromParts[i] !== toParts[i]) {
            samePartsLength = i;
            break;
          }
        }
        var outputParts = [];
        for (var i = samePartsLength; i < fromParts.length; i++) {
          outputParts.push("..");
        }
        outputParts = outputParts.concat(toParts.slice(samePartsLength));
        return outputParts.join("/");
      },
    };
    var FS_stdin_getChar_buffer = [];
    function intArrayFromString(stringy, dontAddNull, length) {
      var len = length > 0 ? length : lengthBytesUTF8(stringy) + 1;
      var u8array = new Array(len);
      var numBytesWritten = stringToUTF8Array(
        stringy,
        u8array,
        0,
        u8array.length,
      );
      if (dontAddNull) u8array.length = numBytesWritten;
      return u8array;
    }
    var FS_stdin_getChar = () => {
      if (!FS_stdin_getChar_buffer.length) {
        var result = null;
        if (ENVIRONMENT_IS_NODE) {
          var BUFSIZE = 256;
          var buf = Buffer.alloc(BUFSIZE);
          var bytesRead = 0;
          var fd = process.stdin.fd;
          try {
            bytesRead = fs.readSync(fd, buf, 0, BUFSIZE, -1);
          } catch (e) {
            if (e.toString().includes("EOF")) bytesRead = 0;
            else throw e;
          }
          if (bytesRead > 0) {
            result = buf.slice(0, bytesRead).toString("utf-8");
          } else {
            result = null;
          }
        } else if (
          typeof window != "undefined" &&
          typeof window.prompt == "function"
        ) {
          result = window.prompt("Input: ");
          if (result !== null) {
            result += "\n";
          }
        } else if (typeof readline == "function") {
          result = readline();
          if (result !== null) {
            result += "\n";
          }
        }
        if (!result) {
          return null;
        }
        FS_stdin_getChar_buffer = intArrayFromString(result, true);
      }
      return FS_stdin_getChar_buffer.shift();
    };
    var TTY = {
      ttys: [],
      init: function () {},
      shutdown: function () {},
      register: function (dev, ops) {
        TTY.ttys[dev] = { input: [], output: [], ops: ops };
        FS.registerDevice(dev, TTY.stream_ops);
      },
      stream_ops: {
        open: function (stream) {
          var tty = TTY.ttys[stream.node.rdev];
          if (!tty) {
            throw new FS.ErrnoError(43);
          }
          stream.tty = tty;
          stream.seekable = false;
        },
        close: function (stream) {
          stream.tty.ops.fsync(stream.tty);
        },
        fsync: function (stream) {
          stream.tty.ops.fsync(stream.tty);
        },
        read: function (stream, buffer, offset, length, pos) {
          if (!stream.tty || !stream.tty.ops.get_char) {
            throw new FS.ErrnoError(60);
          }
          var bytesRead = 0;
          for (var i = 0; i < length; i++) {
            var result;
            try {
              result = stream.tty.ops.get_char(stream.tty);
            } catch (e) {
              throw new FS.ErrnoError(29);
            }
            if (result === undefined && bytesRead === 0) {
              throw new FS.ErrnoError(6);
            }
            if (result === null || result === undefined) break;
            bytesRead++;
            buffer[offset + i] = result;
          }
          if (bytesRead) {
            stream.node.timestamp = Date.now();
          }
          return bytesRead;
        },
        write: function (stream, buffer, offset, length, pos) {
          if (!stream.tty || !stream.tty.ops.put_char) {
            throw new FS.ErrnoError(60);
          }
          try {
            for (var i = 0; i < length; i++) {
              stream.tty.ops.put_char(stream.tty, buffer[offset + i]);
            }
          } catch (e) {
            throw new FS.ErrnoError(29);
          }
          if (length) {
            stream.node.timestamp = Date.now();
          }
          return i;
        },
      },
      default_tty_ops: {
        get_char: function (tty) {
          return FS_stdin_getChar();
        },
        put_char: function (tty, val) {
          if (val === null || val === 10) {
            out(UTF8ArrayToString(tty.output, 0));
            tty.output = [];
          } else {
            if (val != 0) tty.output.push(val);
          }
        },
        fsync: function (tty) {
          if (tty.output && tty.output.length > 0) {
            out(UTF8ArrayToString(tty.output, 0));
            tty.output = [];
          }
        },
        ioctl_tcgets: function (tty) {
          return {
            c_iflag: 25856,
            c_oflag: 5,
            c_cflag: 191,
            c_lflag: 35387,
            c_cc: [
              3, 28, 127, 21, 4, 0, 1, 0, 17, 19, 26, 0, 18, 15, 23, 22, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ],
          };
        },
        ioctl_tcsets: function (tty, optional_actions, data) {
          return 0;
        },
        ioctl_tiocgwinsz: function (tty) {
          return [24, 80];
        },
      },
      default_tty1_ops: {
        put_char: function (tty, val) {
          if (val === null || val === 10) {
            err(UTF8ArrayToString(tty.output, 0));
            tty.output = [];
          } else {
            if (val != 0) tty.output.push(val);
          }
        },
        fsync: function (tty) {
          if (tty.output && tty.output.length > 0) {
            err(UTF8ArrayToString(tty.output, 0));
            tty.output = [];
          }
        },
      },
    };
    var zeroMemory = (address, size) => {
      HEAPU8.fill(0, address, address + size);
      return address;
    };
    var alignMemory = (size, alignment) => {
      assert(alignment, "alignment argument is required");
      return Math.ceil(size / alignment) * alignment;
    };
    var mmapAlloc = (size) => {
      size = alignMemory(size, 65536);
      var ptr = _emscripten_builtin_memalign(65536, size);
      if (!ptr) return 0;
      return zeroMemory(ptr, size);
    };
    var MEMFS = {
      ops_table: null,
      mount(mount) {
        return MEMFS.createNode(null, "/", 16384 | 511, 0);
      },
      createNode(parent, name, mode, dev) {
        if (FS.isBlkdev(mode) || FS.isFIFO(mode)) {
          throw new FS.ErrnoError(63);
        }
        if (!MEMFS.ops_table) {
          MEMFS.ops_table = {
            dir: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
                lookup: MEMFS.node_ops.lookup,
                mknod: MEMFS.node_ops.mknod,
                rename: MEMFS.node_ops.rename,
                unlink: MEMFS.node_ops.unlink,
                rmdir: MEMFS.node_ops.rmdir,
                readdir: MEMFS.node_ops.readdir,
                symlink: MEMFS.node_ops.symlink,
              },
              stream: { llseek: MEMFS.stream_ops.llseek },
            },
            file: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
              },
              stream: {
                llseek: MEMFS.stream_ops.llseek,
                read: MEMFS.stream_ops.read,
                write: MEMFS.stream_ops.write,
                allocate: MEMFS.stream_ops.allocate,
                mmap: MEMFS.stream_ops.mmap,
                msync: MEMFS.stream_ops.msync,
              },
            },
            link: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
                readlink: MEMFS.node_ops.readlink,
              },
              stream: {},
            },
            chrdev: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
              },
              stream: FS.chrdev_stream_ops,
            },
          };
        }
        var node = FS.createNode(parent, name, mode, dev);
        if (FS.isDir(node.mode)) {
          node.node_ops = MEMFS.ops_table.dir.node;
          node.stream_ops = MEMFS.ops_table.dir.stream;
          node.contents = {};
        } else if (FS.isFile(node.mode)) {
          node.node_ops = MEMFS.ops_table.file.node;
          node.stream_ops = MEMFS.ops_table.file.stream;
          node.usedBytes = 0;
          node.contents = null;
        } else if (FS.isLink(node.mode)) {
          node.node_ops = MEMFS.ops_table.link.node;
          node.stream_ops = MEMFS.ops_table.link.stream;
        } else if (FS.isChrdev(node.mode)) {
          node.node_ops = MEMFS.ops_table.chrdev.node;
          node.stream_ops = MEMFS.ops_table.chrdev.stream;
        }
        node.timestamp = Date.now();
        if (parent) {
          parent.contents[name] = node;
          parent.timestamp = node.timestamp;
        }
        return node;
      },
      getFileDataAsTypedArray(node) {
        if (!node.contents) return new Uint8Array(0);
        if (node.contents.subarray)
          return node.contents.subarray(0, node.usedBytes);
        return new Uint8Array(node.contents);
      },
      expandFileStorage(node, newCapacity) {
        var prevCapacity = node.contents ? node.contents.length : 0;
        if (prevCapacity >= newCapacity) return;
        var CAPACITY_DOUBLING_MAX = 1024 * 1024;
        newCapacity = Math.max(
          newCapacity,
          (prevCapacity *
            (prevCapacity < CAPACITY_DOUBLING_MAX ? 2 : 1.125)) >>>
            0,
        );
        if (prevCapacity != 0) newCapacity = Math.max(newCapacity, 256);
        var oldContents = node.contents;
        node.contents = new Uint8Array(newCapacity);
        if (node.usedBytes > 0)
          node.contents.set(oldContents.subarray(0, node.usedBytes), 0);
      },
      resizeFileStorage(node, newSize) {
        if (node.usedBytes == newSize) return;
        if (newSize == 0) {
          node.contents = null;
          node.usedBytes = 0;
        } else {
          var oldContents = node.contents;
          node.contents = new Uint8Array(newSize);
          if (oldContents) {
            node.contents.set(
              oldContents.subarray(0, Math.min(newSize, node.usedBytes)),
            );
          }
          node.usedBytes = newSize;
        }
      },
      node_ops: {
        getattr(node) {
          var attr = {};
          attr.dev = FS.isChrdev(node.mode) ? node.id : 1;
          attr.ino = node.id;
          attr.mode = node.mode;
          attr.nlink = 1;
          attr.uid = 0;
          attr.gid = 0;
          attr.rdev = node.rdev;
          if (FS.isDir(node.mode)) {
            attr.size = 4096;
          } else if (FS.isFile(node.mode)) {
            attr.size = node.usedBytes;
          } else if (FS.isLink(node.mode)) {
            attr.size = node.link.length;
          } else {
            attr.size = 0;
          }
          attr.atime = new Date(node.timestamp);
          attr.mtime = new Date(node.timestamp);
          attr.ctime = new Date(node.timestamp);
          attr.blksize = 4096;
          attr.blocks = Math.ceil(attr.size / attr.blksize);
          return attr;
        },
        setattr(node, attr) {
          if (attr.mode !== undefined) {
            node.mode = attr.mode;
          }
          if (attr.timestamp !== undefined) {
            node.timestamp = attr.timestamp;
          }
          if (attr.size !== undefined) {
            MEMFS.resizeFileStorage(node, attr.size);
          }
        },
        lookup(parent, name) {
          throw FS.genericErrors[44];
        },
        mknod(parent, name, mode, dev) {
          return MEMFS.createNode(parent, name, mode, dev);
        },
        rename(old_node, new_dir, new_name) {
          if (FS.isDir(old_node.mode)) {
            var new_node;
            try {
              new_node = FS.lookupNode(new_dir, new_name);
            } catch (e) {}
            if (new_node) {
              for (var i in new_node.contents) {
                throw new FS.ErrnoError(55);
              }
            }
          }
          delete old_node.parent.contents[old_node.name];
          old_node.parent.timestamp = Date.now();
          old_node.name = new_name;
          new_dir.contents[new_name] = old_node;
          new_dir.timestamp = old_node.parent.timestamp;
          old_node.parent = new_dir;
        },
        unlink(parent, name) {
          delete parent.contents[name];
          parent.timestamp = Date.now();
        },
        rmdir(parent, name) {
          var node = FS.lookupNode(parent, name);
          for (var i in node.contents) {
            throw new FS.ErrnoError(55);
          }
          delete parent.contents[name];
          parent.timestamp = Date.now();
        },
        readdir(node) {
          var entries = [".", ".."];
          for (var key in node.contents) {
            if (!node.contents.hasOwnProperty(key)) {
              continue;
            }
            entries.push(key);
          }
          return entries;
        },
        symlink(parent, newname, oldpath) {
          var node = MEMFS.createNode(parent, newname, 511 | 40960, 0);
          node.link = oldpath;
          return node;
        },
        readlink(node) {
          if (!FS.isLink(node.mode)) {
            throw new FS.ErrnoError(28);
          }
          return node.link;
        },
      },
      stream_ops: {
        read(stream, buffer, offset, length, position) {
          var contents = stream.node.contents;
          if (position >= stream.node.usedBytes) return 0;
          var size = Math.min(stream.node.usedBytes - position, length);
          assert(size >= 0);
          if (size > 8 && contents.subarray) {
            buffer.set(contents.subarray(position, position + size), offset);
          } else {
            for (var i = 0; i < size; i++)
              buffer[offset + i] = contents[position + i];
          }
          return size;
        },
        write(stream, buffer, offset, length, position, canOwn) {
          assert(!(buffer instanceof ArrayBuffer));
          if (buffer.buffer === HEAP8.buffer) {
            canOwn = false;
          }
          if (!length) return 0;
          var node = stream.node;
          node.timestamp = Date.now();
          if (buffer.subarray && (!node.contents || node.contents.subarray)) {
            if (canOwn) {
              assert(
                position === 0,
                "canOwn must imply no weird position inside the file",
              );
              node.contents = buffer.subarray(offset, offset + length);
              node.usedBytes = length;
              return length;
            } else if (node.usedBytes === 0 && position === 0) {
              node.contents = buffer.slice(offset, offset + length);
              node.usedBytes = length;
              return length;
            } else if (position + length <= node.usedBytes) {
              node.contents.set(
                buffer.subarray(offset, offset + length),
                position,
              );
              return length;
            }
          }
          MEMFS.expandFileStorage(node, position + length);
          if (node.contents.subarray && buffer.subarray) {
            node.contents.set(
              buffer.subarray(offset, offset + length),
              position,
            );
          } else {
            for (var i = 0; i < length; i++) {
              node.contents[position + i] = buffer[offset + i];
            }
          }
          node.usedBytes = Math.max(node.usedBytes, position + length);
          return length;
        },
        llseek(stream, offset, whence) {
          var position = offset;
          if (whence === 1) {
            position += stream.position;
          } else if (whence === 2) {
            if (FS.isFile(stream.node.mode)) {
              position += stream.node.usedBytes;
            }
          }
          if (position < 0) {
            throw new FS.ErrnoError(28);
          }
          return position;
        },
        allocate(stream, offset, length) {
          MEMFS.expandFileStorage(stream.node, offset + length);
          stream.node.usedBytes = Math.max(
            stream.node.usedBytes,
            offset + length,
          );
        },
        mmap(stream, length, position, prot, flags) {
          if (!FS.isFile(stream.node.mode)) {
            throw new FS.ErrnoError(43);
          }
          var ptr;
          var allocated;
          var contents = stream.node.contents;
          if (!(flags & 2) && contents.buffer === HEAP8.buffer) {
            allocated = false;
            ptr = contents.byteOffset;
          } else {
            if (position > 0 || position + length < contents.length) {
              if (contents.subarray) {
                contents = contents.subarray(position, position + length);
              } else {
                contents = Array.prototype.slice.call(
                  contents,
                  position,
                  position + length,
                );
              }
            }
            allocated = true;
            ptr = mmapAlloc(length);
            if (!ptr) {
              throw new FS.ErrnoError(48);
            }
            HEAP8.set(contents, ptr);
          }
          return { ptr: ptr, allocated: allocated };
        },
        msync(stream, buffer, offset, length, mmapFlags) {
          MEMFS.stream_ops.write(stream, buffer, 0, length, offset, false);
          return 0;
        },
      },
    };
    var asyncLoad = (url, onload, onerror, noRunDep) => {
      var dep = !noRunDep ? getUniqueRunDependency(`al ${url}`) : "";
      readAsync(
        url,
        (arrayBuffer) => {
          assert(
            arrayBuffer,
            `Loading data file "${url}" failed (no arrayBuffer).`,
          );
          onload(new Uint8Array(arrayBuffer));
          if (dep) removeRunDependency(dep);
        },
        (event) => {
          if (onerror) {
            onerror();
          } else {
            throw `Loading data file "${url}" failed.`;
          }
        },
      );
      if (dep) addRunDependency(dep);
    };
    var preloadPlugins = Module["preloadPlugins"] || [];
    function FS_handledByPreloadPlugin(byteArray, fullname, finish, onerror) {
      if (typeof Browser != "undefined") Browser.init();
      var handled = false;
      preloadPlugins.forEach(function (plugin) {
        if (handled) return;
        if (plugin["canHandle"](fullname)) {
          plugin["handle"](byteArray, fullname, finish, onerror);
          handled = true;
        }
      });
      return handled;
    }
    function FS_createPreloadedFile(
      parent,
      name,
      url,
      canRead,
      canWrite,
      onload,
      onerror,
      dontCreateFile,
      canOwn,
      preFinish,
    ) {
      var fullname = name ? PATH_FS.resolve(PATH.join2(parent, name)) : parent;
      var dep = getUniqueRunDependency(`cp ${fullname}`);
      function processData(byteArray) {
        function finish(byteArray) {
          if (preFinish) preFinish();
          if (!dontCreateFile) {
            FS.createDataFile(
              parent,
              name,
              byteArray,
              canRead,
              canWrite,
              canOwn,
            );
          }
          if (onload) onload();
          removeRunDependency(dep);
        }
        if (
          FS_handledByPreloadPlugin(byteArray, fullname, finish, () => {
            if (onerror) onerror();
            removeRunDependency(dep);
          })
        ) {
          return;
        }
        finish(byteArray);
      }
      addRunDependency(dep);
      if (typeof url == "string") {
        asyncLoad(url, (byteArray) => processData(byteArray), onerror);
      } else {
        processData(url);
      }
    }
    function FS_modeStringToFlags(str) {
      var flagModes = {
        r: 0,
        "r+": 2,
        w: 512 | 64 | 1,
        "w+": 512 | 64 | 2,
        a: 1024 | 64 | 1,
        "a+": 1024 | 64 | 2,
      };
      var flags = flagModes[str];
      if (typeof flags == "undefined") {
        throw new Error(`Unknown file open mode: ${str}`);
      }
      return flags;
    }
    function FS_getMode(canRead, canWrite) {
      var mode = 0;
      if (canRead) mode |= 292 | 73;
      if (canWrite) mode |= 146;
      return mode;
    }
    var ERRNO_CODES = {};
    var NODEFS = {
      isWindows: false,
      staticInit: () => {
        NODEFS.isWindows = !!process.platform.match(/^win/);
        var flags = process.binding("constants");
        if (flags["fs"]) {
          flags = flags["fs"];
        }
        NODEFS.flagsForNodeMap = {
          1024: flags["O_APPEND"],
          64: flags["O_CREAT"],
          128: flags["O_EXCL"],
          256: flags["O_NOCTTY"],
          0: flags["O_RDONLY"],
          2: flags["O_RDWR"],
          4096: flags["O_SYNC"],
          512: flags["O_TRUNC"],
          1: flags["O_WRONLY"],
          131072: flags["O_NOFOLLOW"],
        };
        assert(NODEFS.flagsForNodeMap["0"] === 0);
      },
      convertNodeCode: (e) => {
        var code = e.code;
        assert(
          code in ERRNO_CODES,
          `unexpected node error code: ${code} (${e})`,
        );
        return ERRNO_CODES[code];
      },
      mount: (mount) => {
        assert(ENVIRONMENT_IS_NODE);
        return NODEFS.createNode(null, "/", NODEFS.getMode(mount.opts.root), 0);
      },
      createNode: (parent, name, mode, dev) => {
        if (!FS.isDir(mode) && !FS.isFile(mode) && !FS.isLink(mode)) {
          throw new FS.ErrnoError(28);
        }
        var node = FS.createNode(parent, name, mode);
        node.node_ops = NODEFS.node_ops;
        node.stream_ops = NODEFS.stream_ops;
        return node;
      },
      getMode: (path) => {
        var stat;
        try {
          stat = fs.lstatSync(path);
          if (NODEFS.isWindows) {
            stat.mode = stat.mode | ((stat.mode & 292) >> 2);
          }
        } catch (e) {
          if (!e.code) throw e;
          throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
        }
        return stat.mode;
      },
      realPath: (node) => {
        var parts = [];
        while (node.parent !== node) {
          parts.push(node.name);
          node = node.parent;
        }
        parts.push(node.mount.opts.root);
        parts.reverse();
        return PATH.join.apply(null, parts);
      },
      flagsForNode: (flags) => {
        flags &= ~2097152;
        flags &= ~2048;
        flags &= ~32768;
        flags &= ~524288;
        flags &= ~65536;
        var newFlags = 0;
        for (var k in NODEFS.flagsForNodeMap) {
          if (flags & k) {
            newFlags |= NODEFS.flagsForNodeMap[k];
            flags ^= k;
          }
        }
        if (flags) {
          throw new FS.ErrnoError(28);
        }
        return newFlags;
      },
      node_ops: {
        getattr: (node) => {
          var path = NODEFS.realPath(node);
          var stat;
          try {
            stat = fs.lstatSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
          if (NODEFS.isWindows && !stat.blksize) {
            stat.blksize = 4096;
          }
          if (NODEFS.isWindows && !stat.blocks) {
            stat.blocks = ((stat.size + stat.blksize - 1) / stat.blksize) | 0;
          }
          return {
            dev: stat.dev,
            ino: stat.ino,
            mode: stat.mode,
            nlink: stat.nlink,
            uid: stat.uid,
            gid: stat.gid,
            rdev: stat.rdev,
            size: stat.size,
            atime: stat.atime,
            mtime: stat.mtime,
            ctime: stat.ctime,
            blksize: stat.blksize,
            blocks: stat.blocks,
          };
        },
        setattr: (node, attr) => {
          var path = NODEFS.realPath(node);
          try {
            if (attr.mode !== undefined) {
              fs.chmodSync(path, attr.mode);
              node.mode = attr.mode;
            }
            if (attr.timestamp !== undefined) {
              var date = new Date(attr.timestamp);
              fs.utimesSync(path, date, date);
            }
            if (attr.size !== undefined) {
              fs.truncateSync(path, attr.size);
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        lookup: (parent, name) => {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          var mode = NODEFS.getMode(path);
          return NODEFS.createNode(parent, name, mode);
        },
        mknod: (parent, name, mode, dev) => {
          var node = NODEFS.createNode(parent, name, mode, dev);
          var path = NODEFS.realPath(node);
          try {
            if (FS.isDir(node.mode)) {
              fs.mkdirSync(path, node.mode);
            } else {
              fs.writeFileSync(path, "", { mode: node.mode });
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
          return node;
        },
        rename: (oldNode, newDir, newName) => {
          var oldPath = NODEFS.realPath(oldNode);
          var newPath = PATH.join2(NODEFS.realPath(newDir), newName);
          try {
            fs.renameSync(oldPath, newPath);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
          oldNode.name = newName;
        },
        unlink: (parent, name) => {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          try {
            fs.unlinkSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        rmdir: (parent, name) => {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          try {
            fs.rmdirSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        readdir: (node) => {
          var path = NODEFS.realPath(node);
          try {
            return fs.readdirSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        symlink: (parent, newName, oldPath) => {
          var newPath = PATH.join2(NODEFS.realPath(parent), newName);
          try {
            fs.symlinkSync(oldPath, newPath);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        readlink: (node) => {
          var path = NODEFS.realPath(node);
          try {
            path = fs.readlinkSync(path);
            path = nodePath.relative(
              nodePath.resolve(node.mount.opts.root),
              path,
            );
            return path;
          } catch (e) {
            if (!e.code) throw e;
            if (e.code === "UNKNOWN") throw new FS.ErrnoError(28);
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
      },
      stream_ops: {
        open: (stream) => {
          var path = NODEFS.realPath(stream.node);
          try {
            if (FS.isFile(stream.node.mode)) {
              stream.nfd = fs.openSync(path, NODEFS.flagsForNode(stream.flags));
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        close: (stream) => {
          try {
            if (FS.isFile(stream.node.mode) && stream.nfd) {
              fs.closeSync(stream.nfd);
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        read: (stream, buffer, offset, length, position) => {
          if (length === 0) return 0;
          try {
            return fs.readSync(
              stream.nfd,
              Buffer.from(buffer.buffer),
              offset,
              length,
              position,
            );
          } catch (e) {
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        write: (stream, buffer, offset, length, position) => {
          try {
            return fs.writeSync(
              stream.nfd,
              Buffer.from(buffer.buffer),
              offset,
              length,
              position,
            );
          } catch (e) {
            throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
          }
        },
        llseek: (stream, offset, whence) => {
          var position = offset;
          if (whence === 1) {
            position += stream.position;
          } else if (whence === 2) {
            if (FS.isFile(stream.node.mode)) {
              try {
                var stat = fs.fstatSync(stream.nfd);
                position += stat.size;
              } catch (e) {
                throw new FS.ErrnoError(NODEFS.convertNodeCode(e));
              }
            }
          }
          if (position < 0) {
            throw new FS.ErrnoError(28);
          }
          return position;
        },
        mmap: (stream, length, position, prot, flags) => {
          if (!FS.isFile(stream.node.mode)) {
            throw new FS.ErrnoError(43);
          }
          var ptr = mmapAlloc(length);
          NODEFS.stream_ops.read(stream, HEAP8, ptr, length, position);
          return { ptr: ptr, allocated: true };
        },
        msync: (stream, buffer, offset, length, mmapFlags) => {
          NODEFS.stream_ops.write(stream, buffer, 0, length, offset, false);
          return 0;
        },
      },
    };
    var ERRNO_MESSAGES = {
      0: "Success",
      1: "Arg list too long",
      2: "Permission denied",
      3: "Address already in use",
      4: "Address not available",
      5: "Address family not supported by protocol family",
      6: "No more processes",
      7: "Socket already connected",
      8: "Bad file number",
      9: "Trying to read unreadable message",
      10: "Mount device busy",
      11: "Operation canceled",
      12: "No children",
      13: "Connection aborted",
      14: "Connection refused",
      15: "Connection reset by peer",
      16: "File locking deadlock error",
      17: "Destination address required",
      18: "Math arg out of domain of func",
      19: "Quota exceeded",
      20: "File exists",
      21: "Bad address",
      22: "File too large",
      23: "Host is unreachable",
      24: "Identifier removed",
      25: "Illegal byte sequence",
      26: "Connection already in progress",
      27: "Interrupted system call",
      28: "Invalid argument",
      29: "I/O error",
      30: "Socket is already connected",
      31: "Is a directory",
      32: "Too many symbolic links",
      33: "Too many open files",
      34: "Too many links",
      35: "Message too long",
      36: "Multihop attempted",
      37: "File or path name too long",
      38: "Network interface is not configured",
      39: "Connection reset by network",
      40: "Network is unreachable",
      41: "Too many open files in system",
      42: "No buffer space available",
      43: "No such device",
      44: "No such file or directory",
      45: "Exec format error",
      46: "No record locks available",
      47: "The link has been severed",
      48: "Not enough core",
      49: "No message of desired type",
      50: "Protocol not available",
      51: "No space left on device",
      52: "Function not implemented",
      53: "Socket is not connected",
      54: "Not a directory",
      55: "Directory not empty",
      56: "State not recoverable",
      57: "Socket operation on non-socket",
      59: "Not a typewriter",
      60: "No such device or address",
      61: "Value too large for defined data type",
      62: "Previous owner died",
      63: "Not super-user",
      64: "Broken pipe",
      65: "Protocol error",
      66: "Unknown protocol",
      67: "Protocol wrong type for socket",
      68: "Math result not representable",
      69: "Read only file system",
      70: "Illegal seek",
      71: "No such process",
      72: "Stale file handle",
      73: "Connection timed out",
      74: "Text file busy",
      75: "Cross-device link",
      100: "Device not a stream",
      101: "Bad font file fmt",
      102: "Invalid slot",
      103: "Invalid request code",
      104: "No anode",
      105: "Block device required",
      106: "Channel number out of range",
      107: "Level 3 halted",
      108: "Level 3 reset",
      109: "Link number out of range",
      110: "Protocol driver not attached",
      111: "No CSI structure available",
      112: "Level 2 halted",
      113: "Invalid exchange",
      114: "Invalid request descriptor",
      115: "Exchange full",
      116: "No data (for no delay io)",
      117: "Timer expired",
      118: "Out of streams resources",
      119: "Machine is not on the network",
      120: "Package not installed",
      121: "The object is remote",
      122: "Advertise error",
      123: "Srmount error",
      124: "Communication error on send",
      125: "Cross mount point (not really error)",
      126: "Given log. name not unique",
      127: "f.d. invalid for this operation",
      128: "Remote address changed",
      129: "Can   access a needed shared lib",
      130: "Accessing a corrupted shared lib",
      131: ".lib section in a.out corrupted",
      132: "Attempting to link in too many libs",
      133: "Attempting to exec a shared library",
      135: "Streams pipe error",
      136: "Too many users",
      137: "Socket type not supported",
      138: "Not supported",
      139: "Protocol family not supported",
      140: "Can't send after socket shutdown",
      141: "Too many references",
      142: "Host is down",
      148: "No medium (in tape drive)",
      156: "Level 2 not synchronized",
    };
    var FS = {
      root: null,
      mounts: [],
      devices: {},
      streams: [],
      nextInode: 1,
      nameTable: null,
      currentPath: "/",
      initialized: false,
      ignorePermissions: true,
      ErrnoError: null,
      genericErrors: {},
      filesystems: null,
      syncFSRequests: 0,
      lookupPath: (path, opts = {}) => {
        path = PATH_FS.resolve(path);
        if (!path) return { path: "", node: null };
        var defaults = { follow_mount: true, recurse_count: 0 };
        opts = Object.assign(defaults, opts);
        if (opts.recurse_count > 8) {
          throw new FS.ErrnoError(32);
        }
        var parts = path.split("/").filter((p) => !!p);
        var current = FS.root;
        var current_path = "/";
        for (var i = 0; i < parts.length; i++) {
          var islast = i === parts.length - 1;
          if (islast && opts.parent) {
            break;
          }
          current = FS.lookupNode(current, parts[i]);
          current_path = PATH.join2(current_path, parts[i]);
          if (FS.isMountpoint(current)) {
            if (!islast || (islast && opts.follow_mount)) {
              current = current.mounted.root;
            }
          }
          if (!islast || opts.follow) {
            var count = 0;
            while (FS.isLink(current.mode)) {
              var link = FS.readlink(current_path);
              current_path = PATH_FS.resolve(PATH.dirname(current_path), link);
              var lookup = FS.lookupPath(current_path, {
                recurse_count: opts.recurse_count + 1,
              });
              current = lookup.node;
              if (count++ > 40) {
                throw new FS.ErrnoError(32);
              }
            }
          }
        }
        return { path: current_path, node: current };
      },
      getPath: (node) => {
        var path;
        while (true) {
          if (FS.isRoot(node)) {
            var mount = node.mount.mountpoint;
            if (!path) return mount;
            return mount[mount.length - 1] !== "/"
              ? `${mount}/${path}`
              : mount + path;
          }
          path = path ? `${node.name}/${path}` : node.name;
          node = node.parent;
        }
      },
      hashName: (parentid, name) => {
        var hash = 0;
        for (var i = 0; i < name.length; i++) {
          hash = ((hash << 5) - hash + name.charCodeAt(i)) | 0;
        }
        return ((parentid + hash) >>> 0) % FS.nameTable.length;
      },
      hashAddNode: (node) => {
        var hash = FS.hashName(node.parent.id, node.name);
        node.name_next = FS.nameTable[hash];
        FS.nameTable[hash] = node;
      },
      hashRemoveNode: (node) => {
        var hash = FS.hashName(node.parent.id, node.name);
        if (FS.nameTable[hash] === node) {
          FS.nameTable[hash] = node.name_next;
        } else {
          var current = FS.nameTable[hash];
          while (current) {
            if (current.name_next === node) {
              current.name_next = node.name_next;
              break;
            }
            current = current.name_next;
          }
        }
      },
      lookupNode: (parent, name) => {
        var errCode = FS.mayLookup(parent);
        if (errCode) {
          throw new FS.ErrnoError(errCode, parent);
        }
        var hash = FS.hashName(parent.id, name);
        for (var node = FS.nameTable[hash]; node; node = node.name_next) {
          var nodeName = node.name;
          if (node.parent.id === parent.id && nodeName === name) {
            return node;
          }
        }
        return FS.lookup(parent, name);
      },
      createNode: (parent, name, mode, rdev) => {
        assert(typeof parent == "object");
        var node = new FS.FSNode(parent, name, mode, rdev);
        FS.hashAddNode(node);
        return node;
      },
      destroyNode: (node) => {
        FS.hashRemoveNode(node);
      },
      isRoot: (node) => node === node.parent,
      isMountpoint: (node) => !!node.mounted,
      isFile: (mode) => (mode & 61440) === 32768,
      isDir: (mode) => (mode & 61440) === 16384,
      isLink: (mode) => (mode & 61440) === 40960,
      isChrdev: (mode) => (mode & 61440) === 8192,
      isBlkdev: (mode) => (mode & 61440) === 24576,
      isFIFO: (mode) => (mode & 61440) === 4096,
      isSocket: (mode) => (mode & 49152) === 49152,
      flagsToPermissionString: (flag) => {
        var perms = ["r", "w", "rw"][flag & 3];
        if (flag & 512) {
          perms += "w";
        }
        return perms;
      },
      nodePermissions: (node, perms) => {
        if (FS.ignorePermissions) {
          return 0;
        }
        if (perms.includes("r") && !(node.mode & 292)) {
          return 2;
        } else if (perms.includes("w") && !(node.mode & 146)) {
          return 2;
        } else if (perms.includes("x") && !(node.mode & 73)) {
          return 2;
        }
        return 0;
      },
      mayLookup: (dir) => {
        var errCode = FS.nodePermissions(dir, "x");
        if (errCode) return errCode;
        if (!dir.node_ops.lookup) return 2;
        return 0;
      },
      mayCreate: (dir, name) => {
        try {
          var node = FS.lookupNode(dir, name);
          return 20;
        } catch (e) {}
        return FS.nodePermissions(dir, "wx");
      },
      mayDelete: (dir, name, isdir) => {
        var node;
        try {
          node = FS.lookupNode(dir, name);
        } catch (e) {
          return e.errno;
        }
        var errCode = FS.nodePermissions(dir, "wx");
        if (errCode) {
          return errCode;
        }
        if (isdir) {
          if (!FS.isDir(node.mode)) {
            return 54;
          }
          if (FS.isRoot(node) || FS.getPath(node) === FS.cwd()) {
            return 10;
          }
        } else {
          if (FS.isDir(node.mode)) {
            return 31;
          }
        }
        return 0;
      },
      mayOpen: (node, flags) => {
        if (!node) {
          return 44;
        }
        if (FS.isLink(node.mode)) {
          return 32;
        } else if (FS.isDir(node.mode)) {
          if (FS.flagsToPermissionString(flags) !== "r" || flags & 512) {
            return 31;
          }
        }
        return FS.nodePermissions(node, FS.flagsToPermissionString(flags));
      },
      MAX_OPEN_FDS: 4096,
      nextfd: () => {
        for (var fd = 0; fd <= FS.MAX_OPEN_FDS; fd++) {
          if (!FS.streams[fd]) {
            return fd;
          }
        }
        throw new FS.ErrnoError(33);
      },
      getStreamChecked: (fd) => {
        var stream = FS.getStream(fd);
        if (!stream) {
          throw new FS.ErrnoError(8);
        }
        return stream;
      },
      getStream: (fd) => FS.streams[fd],
      createStream: (stream, fd = -1) => {
        if (!FS.FSStream) {
          FS.FSStream = function () {
            this.shared = {};
          };
          FS.FSStream.prototype = {};
          Object.defineProperties(FS.FSStream.prototype, {
            object: {
              get() {
                return this.node;
              },
              set(val) {
                this.node = val;
              },
            },
            isRead: {
              get() {
                return (this.flags & 2097155) !== 1;
              },
            },
            isWrite: {
              get() {
                return (this.flags & 2097155) !== 0;
              },
            },
            isAppend: {
              get() {
                return this.flags & 1024;
              },
            },
            flags: {
              get() {
                return this.shared.flags;
              },
              set(val) {
                this.shared.flags = val;
              },
            },
            position: {
              get() {
                return this.shared.position;
              },
              set(val) {
                this.shared.position = val;
              },
            },
          });
        }
        stream = Object.assign(new FS.FSStream(), stream);
        if (fd == -1) {
          fd = FS.nextfd();
        }
        stream.fd = fd;
        FS.streams[fd] = stream;
        return stream;
      },
      closeStream: (fd) => {
        FS.streams[fd] = null;
      },
      chrdev_stream_ops: {
        open: (stream) => {
          var device = FS.getDevice(stream.node.rdev);
          stream.stream_ops = device.stream_ops;
          if (stream.stream_ops.open) {
            stream.stream_ops.open(stream);
          }
        },
        llseek: () => {
          throw new FS.ErrnoError(70);
        },
      },
      major: (dev) => dev >> 8,
      minor: (dev) => dev & 255,
      makedev: (ma, mi) => (ma << 8) | mi,
      registerDevice: (dev, ops) => {
        FS.devices[dev] = { stream_ops: ops };
      },
      getDevice: (dev) => FS.devices[dev],
      getMounts: (mount) => {
        var mounts = [];
        var check = [mount];
        while (check.length) {
          var m = check.pop();
          mounts.push(m);
          check.push.apply(check, m.mounts);
        }
        return mounts;
      },
      syncfs: (populate, callback) => {
        if (typeof populate == "function") {
          callback = populate;
          populate = false;
        }
        FS.syncFSRequests++;
        if (FS.syncFSRequests > 1) {
          err(
            `warning: ${FS.syncFSRequests} FS.syncfs operations in flight at once, probably just doing extra work`,
          );
        }
        var mounts = FS.getMounts(FS.root.mount);
        var completed = 0;
        function doCallback(errCode) {
          assert(FS.syncFSRequests > 0);
          FS.syncFSRequests--;
          return callback(errCode);
        }
        function done(errCode) {
          if (errCode) {
            if (!done.errored) {
              done.errored = true;
              return doCallback(errCode);
            }
            return;
          }
          if (++completed >= mounts.length) {
            doCallback(null);
          }
        }
        mounts.forEach((mount) => {
          if (!mount.type.syncfs) {
            return done(null);
          }
          mount.type.syncfs(mount, populate, done);
        });
      },
      mount: (type, opts, mountpoint) => {
        if (typeof type == "string") {
          throw type;
        }
        var root = mountpoint === "/";
        var pseudo = !mountpoint;
        var node;
        if (root && FS.root) {
          throw new FS.ErrnoError(10);
        } else if (!root && !pseudo) {
          var lookup = FS.lookupPath(mountpoint, { follow_mount: false });
          mountpoint = lookup.path;
          node = lookup.node;
          if (FS.isMountpoint(node)) {
            throw new FS.ErrnoError(10);
          }
          if (!FS.isDir(node.mode)) {
            throw new FS.ErrnoError(54);
          }
        }
        var mount = {
          type: type,
          opts: opts,
          mountpoint: mountpoint,
          mounts: [],
        };
        var mountRoot = type.mount(mount);
        mountRoot.mount = mount;
        mount.root = mountRoot;
        if (root) {
          FS.root = mountRoot;
        } else if (node) {
          node.mounted = mount;
          if (node.mount) {
            node.mount.mounts.push(mount);
          }
        }
        return mountRoot;
      },
      unmount: (mountpoint) => {
        var lookup = FS.lookupPath(mountpoint, { follow_mount: false });
        if (!FS.isMountpoint(lookup.node)) {
          throw new FS.ErrnoError(28);
        }
        var node = lookup.node;
        var mount = node.mounted;
        var mounts = FS.getMounts(mount);
        Object.keys(FS.nameTable).forEach((hash) => {
          var current = FS.nameTable[hash];
          while (current) {
            var next = current.name_next;
            if (mounts.includes(current.mount)) {
              FS.destroyNode(current);
            }
            current = next;
          }
        });
        node.mounted = null;
        var idx = node.mount.mounts.indexOf(mount);
        assert(idx !== -1);
        node.mount.mounts.splice(idx, 1);
      },
      lookup: (parent, name) => parent.node_ops.lookup(parent, name),
      mknod: (path, mode, dev) => {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        if (!name || name === "." || name === "..") {
          throw new FS.ErrnoError(28);
        }
        var errCode = FS.mayCreate(parent, name);
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        if (!parent.node_ops.mknod) {
          throw new FS.ErrnoError(63);
        }
        return parent.node_ops.mknod(parent, name, mode, dev);
      },
      create: (path, mode) => {
        mode = mode !== undefined ? mode : 438;
        mode &= 4095;
        mode |= 32768;
        return FS.mknod(path, mode, 0);
      },
      mkdir: (path, mode) => {
        mode = mode !== undefined ? mode : 511;
        mode &= 511 | 512;
        mode |= 16384;
        return FS.mknod(path, mode, 0);
      },
      mkdirTree: (path, mode) => {
        var dirs = path.split("/");
        var d = "";
        for (var i = 0; i < dirs.length; ++i) {
          if (!dirs[i]) continue;
          d += "/" + dirs[i];
          try {
            FS.mkdir(d, mode);
          } catch (e) {
            if (e.errno != 20) throw e;
          }
        }
      },
      mkdev: (path, mode, dev) => {
        if (typeof dev == "undefined") {
          dev = mode;
          mode = 438;
        }
        mode |= 8192;
        return FS.mknod(path, mode, dev);
      },
      symlink: (oldpath, newpath) => {
        if (!PATH_FS.resolve(oldpath)) {
          throw new FS.ErrnoError(44);
        }
        var lookup = FS.lookupPath(newpath, { parent: true });
        var parent = lookup.node;
        if (!parent) {
          throw new FS.ErrnoError(44);
        }
        var newname = PATH.basename(newpath);
        var errCode = FS.mayCreate(parent, newname);
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        if (!parent.node_ops.symlink) {
          throw new FS.ErrnoError(63);
        }
        return parent.node_ops.symlink(parent, newname, oldpath);
      },
      rename: (old_path, new_path) => {
        var old_dirname = PATH.dirname(old_path);
        var new_dirname = PATH.dirname(new_path);
        var old_name = PATH.basename(old_path);
        var new_name = PATH.basename(new_path);
        var lookup, old_dir, new_dir;
        lookup = FS.lookupPath(old_path, { parent: true });
        old_dir = lookup.node;
        lookup = FS.lookupPath(new_path, { parent: true });
        new_dir = lookup.node;
        if (!old_dir || !new_dir) throw new FS.ErrnoError(44);
        if (old_dir.mount !== new_dir.mount) {
          throw new FS.ErrnoError(75);
        }
        var old_node = FS.lookupNode(old_dir, old_name);
        var relative = PATH_FS.relative(old_path, new_dirname);
        if (relative.charAt(0) !== ".") {
          throw new FS.ErrnoError(28);
        }
        relative = PATH_FS.relative(new_path, old_dirname);
        if (relative.charAt(0) !== ".") {
          throw new FS.ErrnoError(55);
        }
        var new_node;
        try {
          new_node = FS.lookupNode(new_dir, new_name);
        } catch (e) {}
        if (old_node === new_node) {
          return;
        }
        var isdir = FS.isDir(old_node.mode);
        var errCode = FS.mayDelete(old_dir, old_name, isdir);
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        errCode = new_node
          ? FS.mayDelete(new_dir, new_name, isdir)
          : FS.mayCreate(new_dir, new_name);
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        if (!old_dir.node_ops.rename) {
          throw new FS.ErrnoError(63);
        }
        if (
          FS.isMountpoint(old_node) ||
          (new_node && FS.isMountpoint(new_node))
        ) {
          throw new FS.ErrnoError(10);
        }
        if (new_dir !== old_dir) {
          errCode = FS.nodePermissions(old_dir, "w");
          if (errCode) {
            throw new FS.ErrnoError(errCode);
          }
        }
        FS.hashRemoveNode(old_node);
        try {
          old_dir.node_ops.rename(old_node, new_dir, new_name);
        } catch (e) {
          throw e;
        } finally {
          FS.hashAddNode(old_node);
        }
      },
      rmdir: (path) => {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        var node = FS.lookupNode(parent, name);
        var errCode = FS.mayDelete(parent, name, true);
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        if (!parent.node_ops.rmdir) {
          throw new FS.ErrnoError(63);
        }
        if (FS.isMountpoint(node)) {
          throw new FS.ErrnoError(10);
        }
        parent.node_ops.rmdir(parent, name);
        FS.destroyNode(node);
      },
      readdir: (path) => {
        var lookup = FS.lookupPath(path, { follow: true });
        var node = lookup.node;
        if (!node.node_ops.readdir) {
          throw new FS.ErrnoError(54);
        }
        return node.node_ops.readdir(node);
      },
      unlink: (path) => {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        if (!parent) {
          throw new FS.ErrnoError(44);
        }
        var name = PATH.basename(path);
        var node = FS.lookupNode(parent, name);
        var errCode = FS.mayDelete(parent, name, false);
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        if (!parent.node_ops.unlink) {
          throw new FS.ErrnoError(63);
        }
        if (FS.isMountpoint(node)) {
          throw new FS.ErrnoError(10);
        }
        parent.node_ops.unlink(parent, name);
        FS.destroyNode(node);
      },
      readlink: (path) => {
        var lookup = FS.lookupPath(path);
        var link = lookup.node;
        if (!link) {
          throw new FS.ErrnoError(44);
        }
        if (!link.node_ops.readlink) {
          throw new FS.ErrnoError(28);
        }
        return PATH_FS.resolve(
          FS.getPath(link.parent),
          link.node_ops.readlink(link),
        );
      },
      stat: (path, dontFollow) => {
        var lookup = FS.lookupPath(path, { follow: !dontFollow });
        var node = lookup.node;
        if (!node) {
          throw new FS.ErrnoError(44);
        }
        if (!node.node_ops.getattr) {
          throw new FS.ErrnoError(63);
        }
        return node.node_ops.getattr(node);
      },
      lstat: (path) => FS.stat(path, true),
      chmod: (path, mode, dontFollow) => {
        var node;
        if (typeof path == "string") {
          var lookup = FS.lookupPath(path, { follow: !dontFollow });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(63);
        }
        node.node_ops.setattr(node, {
          mode: (mode & 4095) | (node.mode & ~4095),
          timestamp: Date.now(),
        });
      },
      lchmod: (path, mode) => {
        FS.chmod(path, mode, true);
      },
      fchmod: (fd, mode) => {
        var stream = FS.getStreamChecked(fd);
        FS.chmod(stream.node, mode);
      },
      chown: (path, uid, gid, dontFollow) => {
        var node;
        if (typeof path == "string") {
          var lookup = FS.lookupPath(path, { follow: !dontFollow });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(63);
        }
        node.node_ops.setattr(node, { timestamp: Date.now() });
      },
      lchown: (path, uid, gid) => {
        FS.chown(path, uid, gid, true);
      },
      fchown: (fd, uid, gid) => {
        var stream = FS.getStreamChecked(fd);
        FS.chown(stream.node, uid, gid);
      },
      truncate: (path, len) => {
        if (len < 0) {
          throw new FS.ErrnoError(28);
        }
        var node;
        if (typeof path == "string") {
          var lookup = FS.lookupPath(path, { follow: true });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(63);
        }
        if (FS.isDir(node.mode)) {
          throw new FS.ErrnoError(31);
        }
        if (!FS.isFile(node.mode)) {
          throw new FS.ErrnoError(28);
        }
        var errCode = FS.nodePermissions(node, "w");
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        node.node_ops.setattr(node, { size: len, timestamp: Date.now() });
      },
      ftruncate: (fd, len) => {
        var stream = FS.getStreamChecked(fd);
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(28);
        }
        FS.truncate(stream.node, len);
      },
      utime: (path, atime, mtime) => {
        var lookup = FS.lookupPath(path, { follow: true });
        var node = lookup.node;
        node.node_ops.setattr(node, { timestamp: Math.max(atime, mtime) });
      },
      open: (path, flags, mode) => {
        if (path === "") {
          throw new FS.ErrnoError(44);
        }
        flags = typeof flags == "string" ? FS_modeStringToFlags(flags) : flags;
        mode = typeof mode == "undefined" ? 438 : mode;
        if (flags & 64) {
          mode = (mode & 4095) | 32768;
        } else {
          mode = 0;
        }
        var node;
        if (typeof path == "object") {
          node = path;
        } else {
          path = PATH.normalize(path);
          try {
            var lookup = FS.lookupPath(path, { follow: !(flags & 131072) });
            node = lookup.node;
          } catch (e) {}
        }
        var created = false;
        if (flags & 64) {
          if (node) {
            if (flags & 128) {
              throw new FS.ErrnoError(20);
            }
          } else {
            node = FS.mknod(path, mode, 0);
            created = true;
          }
        }
        if (!node) {
          throw new FS.ErrnoError(44);
        }
        if (FS.isChrdev(node.mode)) {
          flags &= ~512;
        }
        if (flags & 65536 && !FS.isDir(node.mode)) {
          throw new FS.ErrnoError(54);
        }
        if (!created) {
          var errCode = FS.mayOpen(node, flags);
          if (errCode) {
            throw new FS.ErrnoError(errCode);
          }
        }
        if (flags & 512 && !created) {
          FS.truncate(node, 0);
        }
        flags &= ~(128 | 512 | 131072);
        var stream = FS.createStream({
          node: node,
          path: FS.getPath(node),
          flags: flags,
          seekable: true,
          position: 0,
          stream_ops: node.stream_ops,
          ungotten: [],
          error: false,
        });
        if (stream.stream_ops.open) {
          stream.stream_ops.open(stream);
        }
        if (Module["logReadFiles"] && !(flags & 1)) {
          if (!FS.readFiles) FS.readFiles = {};
          if (!(path in FS.readFiles)) {
            FS.readFiles[path] = 1;
          }
        }
        return stream;
      },
      close: (stream) => {
        if (FS.isClosed(stream)) {
          throw new FS.ErrnoError(8);
        }
        if (stream.getdents) stream.getdents = null;
        try {
          if (stream.stream_ops.close) {
            stream.stream_ops.close(stream);
          }
        } catch (e) {
          throw e;
        } finally {
          FS.closeStream(stream.fd);
        }
        stream.fd = null;
      },
      isClosed: (stream) => stream.fd === null,
      llseek: (stream, offset, whence) => {
        if (FS.isClosed(stream)) {
          throw new FS.ErrnoError(8);
        }
        if (!stream.seekable || !stream.stream_ops.llseek) {
          throw new FS.ErrnoError(70);
        }
        if (whence != 0 && whence != 1 && whence != 2) {
          throw new FS.ErrnoError(28);
        }
        stream.position = stream.stream_ops.llseek(stream, offset, whence);
        stream.ungotten = [];
        return stream.position;
      },
      read: (stream, buffer, offset, length, position) => {
        assert(offset >= 0);
        if (length < 0 || position < 0) {
          throw new FS.ErrnoError(28);
        }
        if (FS.isClosed(stream)) {
          throw new FS.ErrnoError(8);
        }
        if ((stream.flags & 2097155) === 1) {
          throw new FS.ErrnoError(8);
        }
        if (FS.isDir(stream.node.mode)) {
          throw new FS.ErrnoError(31);
        }
        if (!stream.stream_ops.read) {
          throw new FS.ErrnoError(28);
        }
        var seeking = typeof position != "undefined";
        if (!seeking) {
          position = stream.position;
        } else if (!stream.seekable) {
          throw new FS.ErrnoError(70);
        }
        var bytesRead = stream.stream_ops.read(
          stream,
          buffer,
          offset,
          length,
          position,
        );
        if (!seeking) stream.position += bytesRead;
        return bytesRead;
      },
      write: (stream, buffer, offset, length, position, canOwn) => {
        assert(offset >= 0);
        if (length < 0 || position < 0) {
          throw new FS.ErrnoError(28);
        }
        if (FS.isClosed(stream)) {
          throw new FS.ErrnoError(8);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(8);
        }
        if (FS.isDir(stream.node.mode)) {
          throw new FS.ErrnoError(31);
        }
        if (!stream.stream_ops.write) {
          throw new FS.ErrnoError(28);
        }
        if (stream.seekable && stream.flags & 1024) {
          FS.llseek(stream, 0, 2);
        }
        var seeking = typeof position != "undefined";
        if (!seeking) {
          position = stream.position;
        } else if (!stream.seekable) {
          throw new FS.ErrnoError(70);
        }
        var bytesWritten = stream.stream_ops.write(
          stream,
          buffer,
          offset,
          length,
          position,
          canOwn,
        );
        if (!seeking) stream.position += bytesWritten;
        return bytesWritten;
      },
      allocate: (stream, offset, length) => {
        if (FS.isClosed(stream)) {
          throw new FS.ErrnoError(8);
        }
        if (offset < 0 || length <= 0) {
          throw new FS.ErrnoError(28);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(8);
        }
        if (!FS.isFile(stream.node.mode) && !FS.isDir(stream.node.mode)) {
          throw new FS.ErrnoError(43);
        }
        if (!stream.stream_ops.allocate) {
          throw new FS.ErrnoError(138);
        }
        stream.stream_ops.allocate(stream, offset, length);
      },
      mmap: (stream, length, position, prot, flags) => {
        if (
          (prot & 2) !== 0 &&
          (flags & 2) === 0 &&
          (stream.flags & 2097155) !== 2
        ) {
          throw new FS.ErrnoError(2);
        }
        if ((stream.flags & 2097155) === 1) {
          throw new FS.ErrnoError(2);
        }
        if (!stream.stream_ops.mmap) {
          throw new FS.ErrnoError(43);
        }
        return stream.stream_ops.mmap(stream, length, position, prot, flags);
      },
      msync: (stream, buffer, offset, length, mmapFlags) => {
        assert(offset >= 0);
        if (!stream.stream_ops.msync) {
          return 0;
        }
        return stream.stream_ops.msync(
          stream,
          buffer,
          offset,
          length,
          mmapFlags,
        );
      },
      munmap: (stream) => 0,
      ioctl: (stream, cmd, arg) => {
        if (!stream.stream_ops.ioctl) {
          throw new FS.ErrnoError(59);
        }
        return stream.stream_ops.ioctl(stream, cmd, arg);
      },
      readFile: (path, opts = {}) => {
        opts.flags = opts.flags || 0;
        opts.encoding = opts.encoding || "binary";
        if (opts.encoding !== "utf8" && opts.encoding !== "binary") {
          throw new Error(`Invalid encoding type "${opts.encoding}"`);
        }
        var ret;
        var stream = FS.open(path, opts.flags);
        var stat = FS.stat(path);
        var length = stat.size;
        var buf = new Uint8Array(length);
        FS.read(stream, buf, 0, length, 0);
        if (opts.encoding === "utf8") {
          ret = UTF8ArrayToString(buf, 0);
        } else if (opts.encoding === "binary") {
          ret = buf;
        }
        FS.close(stream);
        return ret;
      },
      writeFile: (path, data, opts = {}) => {
        opts.flags = opts.flags || 577;
        var stream = FS.open(path, opts.flags, opts.mode);
        if (typeof data == "string") {
          var buf = new Uint8Array(lengthBytesUTF8(data) + 1);
          var actualNumBytes = stringToUTF8Array(data, buf, 0, buf.length);
          FS.write(stream, buf, 0, actualNumBytes, undefined, opts.canOwn);
        } else if (ArrayBuffer.isView(data)) {
          FS.write(stream, data, 0, data.byteLength, undefined, opts.canOwn);
        } else {
          throw new Error("Unsupported data type");
        }
        FS.close(stream);
      },
      cwd: () => FS.currentPath,
      chdir: (path) => {
        var lookup = FS.lookupPath(path, { follow: true });
        if (lookup.node === null) {
          throw new FS.ErrnoError(44);
        }
        if (!FS.isDir(lookup.node.mode)) {
          throw new FS.ErrnoError(54);
        }
        var errCode = FS.nodePermissions(lookup.node, "x");
        if (errCode) {
          throw new FS.ErrnoError(errCode);
        }
        FS.currentPath = lookup.path;
      },
      createDefaultDirectories: () => {
        FS.mkdir("/tmp");
        FS.mkdir("/home");
        FS.mkdir("/home/web_user");
      },
      createDefaultDevices: () => {
        FS.mkdir("/dev");
        FS.registerDevice(FS.makedev(1, 3), {
          read: () => 0,
          write: (stream, buffer, offset, length, pos) => length,
        });
        FS.mkdev("/dev/null", FS.makedev(1, 3));
        TTY.register(FS.makedev(5, 0), TTY.default_tty_ops);
        TTY.register(FS.makedev(6, 0), TTY.default_tty1_ops);
        FS.mkdev("/dev/tty", FS.makedev(5, 0));
        FS.mkdev("/dev/tty1", FS.makedev(6, 0));
        var randomBuffer = new Uint8Array(1024),
          randomLeft = 0;
        var randomByte = () => {
          if (randomLeft === 0) {
            randomLeft = randomFill(randomBuffer).byteLength;
          }
          return randomBuffer[--randomLeft];
        };
        FS.createDevice("/dev", "random", randomByte);
        FS.createDevice("/dev", "urandom", randomByte);
        FS.mkdir("/dev/shm");
        FS.mkdir("/dev/shm/tmp");
      },
      createSpecialDirectories: () => {
        FS.mkdir("/proc");
        var proc_self = FS.mkdir("/proc/self");
        FS.mkdir("/proc/self/fd");
        FS.mount(
          {
            mount: () => {
              var node = FS.createNode(proc_self, "fd", 16384 | 511, 73);
              node.node_ops = {
                lookup: (parent, name) => {
                  var fd = +name;
                  var stream = FS.getStreamChecked(fd);
                  var ret = {
                    parent: null,
                    mount: { mountpoint: "fake" },
                    node_ops: { readlink: () => stream.path },
                  };
                  ret.parent = ret;
                  return ret;
                },
              };
              return node;
            },
          },
          {},
          "/proc/self/fd",
        );
      },
      createStandardStreams: () => {
        if (Module["stdin"]) {
          FS.createDevice("/dev", "stdin", Module["stdin"]);
        } else {
          FS.symlink("/dev/tty", "/dev/stdin");
        }
        if (Module["stdout"]) {
          FS.createDevice("/dev", "stdout", null, Module["stdout"]);
        } else {
          FS.symlink("/dev/tty", "/dev/stdout");
        }
        if (Module["stderr"]) {
          FS.createDevice("/dev", "stderr", null, Module["stderr"]);
        } else {
          FS.symlink("/dev/tty1", "/dev/stderr");
        }
        var stdin = FS.open("/dev/stdin", 0);
        var stdout = FS.open("/dev/stdout", 1);
        var stderr = FS.open("/dev/stderr", 1);
        assert(stdin.fd === 0, `invalid handle for stdin (${stdin.fd})`);
        assert(stdout.fd === 1, `invalid handle for stdout (${stdout.fd})`);
        assert(stderr.fd === 2, `invalid handle for stderr (${stderr.fd})`);
      },
      ensureErrnoError: () => {
        if (FS.ErrnoError) return;
        FS.ErrnoError = function ErrnoError(errno, node) {
          this.name = "ErrnoError";
          this.node = node;
          this.setErrno = function (errno) {
            this.errno = errno;
            for (var key in ERRNO_CODES) {
              if (ERRNO_CODES[key] === errno) {
                this.code = key;
                break;
              }
            }
          };
          this.setErrno(errno);
          this.message = ERRNO_MESSAGES[errno];
          if (this.stack) {
            Object.defineProperty(this, "stack", {
              value: new Error().stack,
              writable: true,
            });
            this.stack = demangleAll(this.stack);
          }
        };
        FS.ErrnoError.prototype = new Error();
        FS.ErrnoError.prototype.constructor = FS.ErrnoError;
        [44].forEach((code) => {
          FS.genericErrors[code] = new FS.ErrnoError(code);
          FS.genericErrors[code].stack = "<generic error, no stack>";
        });
      },
      staticInit: () => {
        FS.ensureErrnoError();
        FS.nameTable = new Array(4096);
        FS.mount(MEMFS, {}, "/");
        FS.createDefaultDirectories();
        FS.createDefaultDevices();
        FS.createSpecialDirectories();
        FS.filesystems = { MEMFS: MEMFS, NODEFS: NODEFS };
      },
      init: (input, output, error) => {
        assert(
          !FS.init.initialized,
          "FS.init was previously called. If you want to initialize later with custom parameters, remove any earlier calls (note that one is automatically added to the generated code)",
        );
        FS.init.initialized = true;
        FS.ensureErrnoError();
        Module["stdin"] = input || Module["stdin"];
        Module["stdout"] = output || Module["stdout"];
        Module["stderr"] = error || Module["stderr"];
        FS.createStandardStreams();
      },
      quit: () => {
        FS.init.initialized = false;
        _fflush(0);
        for (var i = 0; i < FS.streams.length; i++) {
          var stream = FS.streams[i];
          if (!stream) {
            continue;
          }
          FS.close(stream);
        }
      },
      findObject: (path, dontResolveLastLink) => {
        var ret = FS.analyzePath(path, dontResolveLastLink);
        if (!ret.exists) {
          return null;
        }
        return ret.object;
      },
      analyzePath: (path, dontResolveLastLink) => {
        try {
          var lookup = FS.lookupPath(path, { follow: !dontResolveLastLink });
          path = lookup.path;
        } catch (e) {}
        var ret = {
          isRoot: false,
          exists: false,
          error: 0,
          name: null,
          path: null,
          object: null,
          parentExists: false,
          parentPath: null,
          parentObject: null,
        };
        try {
          var lookup = FS.lookupPath(path, { parent: true });
          ret.parentExists = true;
          ret.parentPath = lookup.path;
          ret.parentObject = lookup.node;
          ret.name = PATH.basename(path);
          lookup = FS.lookupPath(path, { follow: !dontResolveLastLink });
          ret.exists = true;
          ret.path = lookup.path;
          ret.object = lookup.node;
          ret.name = lookup.node.name;
          ret.isRoot = lookup.path === "/";
        } catch (e) {
          ret.error = e.errno;
        }
        return ret;
      },
      createPath: (parent, path, canRead, canWrite) => {
        parent = typeof parent == "string" ? parent : FS.getPath(parent);
        var parts = path.split("/").reverse();
        while (parts.length) {
          var part = parts.pop();
          if (!part) continue;
          var current = PATH.join2(parent, part);
          try {
            FS.mkdir(current);
          } catch (e) {}
          parent = current;
        }
        return current;
      },
      createFile: (parent, name, properties, canRead, canWrite) => {
        var path = PATH.join2(
          typeof parent == "string" ? parent : FS.getPath(parent),
          name,
        );
        var mode = FS_getMode(canRead, canWrite);
        return FS.create(path, mode);
      },
      createDataFile: (parent, name, data, canRead, canWrite, canOwn) => {
        var path = name;
        if (parent) {
          parent = typeof parent == "string" ? parent : FS.getPath(parent);
          path = name ? PATH.join2(parent, name) : parent;
        }
        var mode = FS_getMode(canRead, canWrite);
        var node = FS.create(path, mode);
        if (data) {
          if (typeof data == "string") {
            var arr = new Array(data.length);
            for (var i = 0, len = data.length; i < len; ++i)
              arr[i] = data.charCodeAt(i);
            data = arr;
          }
          FS.chmod(node, mode | 146);
          var stream = FS.open(node, 577);
          FS.write(stream, data, 0, data.length, 0, canOwn);
          FS.close(stream);
          FS.chmod(node, mode);
        }
        return node;
      },
      createDevice: (parent, name, input, output) => {
        var path = PATH.join2(
          typeof parent == "string" ? parent : FS.getPath(parent),
          name,
        );
        var mode = FS_getMode(!!input, !!output);
        if (!FS.createDevice.major) FS.createDevice.major = 64;
        var dev = FS.makedev(FS.createDevice.major++, 0);
        FS.registerDevice(dev, {
          open: (stream) => {
            stream.seekable = false;
          },
          close: (stream) => {
            if (output && output.buffer && output.buffer.length) {
              output(10);
            }
          },
          read: (stream, buffer, offset, length, pos) => {
            var bytesRead = 0;
            for (var i = 0; i < length; i++) {
              var result;
              try {
                result = input();
              } catch (e) {
                throw new FS.ErrnoError(29);
              }
              if (result === undefined && bytesRead === 0) {
                throw new FS.ErrnoError(6);
              }
              if (result === null || result === undefined) break;
              bytesRead++;
              buffer[offset + i] = result;
            }
            if (bytesRead) {
              stream.node.timestamp = Date.now();
            }
            return bytesRead;
          },
          write: (stream, buffer, offset, length, pos) => {
            for (var i = 0; i < length; i++) {
              try {
                output(buffer[offset + i]);
              } catch (e) {
                throw new FS.ErrnoError(29);
              }
            }
            if (length) {
              stream.node.timestamp = Date.now();
            }
            return i;
          },
        });
        return FS.mkdev(path, mode, dev);
      },
      forceLoadFile: (obj) => {
        if (obj.isDevice || obj.isFolder || obj.link || obj.contents)
          return true;
        if (typeof XMLHttpRequest != "undefined") {
          throw new Error(
            "Lazy loading should have been performed (contents set) in createLazyFile, but it was not. Lazy loading only works in web workers. Use --embed-file or --preload-file in emcc on the main thread.",
          );
        } else if (read_) {
          try {
            obj.contents = intArrayFromString(read_(obj.url), true);
            obj.usedBytes = obj.contents.length;
          } catch (e) {
            throw new FS.ErrnoError(29);
          }
        } else {
          throw new Error("Cannot load without read() or XMLHttpRequest.");
        }
      },
      createLazyFile: (parent, name, url, canRead, canWrite) => {
        function LazyUint8Array() {
          this.lengthKnown = false;
          this.chunks = [];
        }
        LazyUint8Array.prototype.get = function LazyUint8Array_get(idx) {
          if (idx > this.length - 1 || idx < 0) {
            return undefined;
          }
          var chunkOffset = idx % this.chunkSize;
          var chunkNum = (idx / this.chunkSize) | 0;
          return this.getter(chunkNum)[chunkOffset];
        };
        LazyUint8Array.prototype.setDataGetter =
          function LazyUint8Array_setDataGetter(getter) {
            this.getter = getter;
          };
        LazyUint8Array.prototype.cacheLength =
          function LazyUint8Array_cacheLength() {
            var xhr = new XMLHttpRequest();
            xhr.open("HEAD", url, false);
            xhr.send(null);
            if (
              !((xhr.status >= 200 && xhr.status < 300) || xhr.status === 304)
            )
              throw new Error(
                "Couldn't load " + url + ". Status: " + xhr.status,
              );
            var datalength = Number(xhr.getResponseHeader("Content-length"));
            var header;
            var hasByteServing =
              (header = xhr.getResponseHeader("Accept-Ranges")) &&
              header === "bytes";
            var usesGzip =
              (header = xhr.getResponseHeader("Content-Encoding")) &&
              header === "gzip";
            var chunkSize = 1024 * 1024;
            if (!hasByteServing) chunkSize = datalength;
            var doXHR = (from, to) => {
              if (from > to)
                throw new Error(
                  "invalid range (" +
                    from +
                    ", " +
                    to +
                    ") or no bytes requested!",
                );
              if (to > datalength - 1)
                throw new Error(
                  "only " + datalength + " bytes available! programmer error!",
                );
              var xhr = new XMLHttpRequest();
              xhr.open("GET", url, false);
              if (datalength !== chunkSize)
                xhr.setRequestHeader("Range", "bytes=" + from + "-" + to);
              xhr.responseType = "arraybuffer";
              if (xhr.overrideMimeType) {
                xhr.overrideMimeType("text/plain; charset=x-user-defined");
              }
              xhr.send(null);
              if (
                !((xhr.status >= 200 && xhr.status < 300) || xhr.status === 304)
              )
                throw new Error(
                  "Couldn't load " + url + ". Status: " + xhr.status,
                );
              if (xhr.response !== undefined) {
                return new Uint8Array(xhr.response || []);
              }
              return intArrayFromString(xhr.responseText || "", true);
            };
            var lazyArray = this;
            lazyArray.setDataGetter((chunkNum) => {
              var start = chunkNum * chunkSize;
              var end = (chunkNum + 1) * chunkSize - 1;
              end = Math.min(end, datalength - 1);
              if (typeof lazyArray.chunks[chunkNum] == "undefined") {
                lazyArray.chunks[chunkNum] = doXHR(start, end);
              }
              if (typeof lazyArray.chunks[chunkNum] == "undefined")
                throw new Error("doXHR failed!");
              return lazyArray.chunks[chunkNum];
            });
            if (usesGzip || !datalength) {
              chunkSize = datalength = 1;
              datalength = this.getter(0).length;
              chunkSize = datalength;
              out(
                "LazyFiles on gzip forces download of the whole file when length is accessed",
              );
            }
            this._length = datalength;
            this._chunkSize = chunkSize;
            this.lengthKnown = true;
          };
        if (typeof XMLHttpRequest != "undefined") {
          if (!ENVIRONMENT_IS_WORKER)
            throw "Cannot do synchronous binary XHRs outside webworkers in modern browsers. Use --embed-file or --preload-file in emcc";
          var lazyArray = new LazyUint8Array();
          Object.defineProperties(lazyArray, {
            length: {
              get: function () {
                if (!this.lengthKnown) {
                  this.cacheLength();
                }
                return this._length;
              },
            },
            chunkSize: {
              get: function () {
                if (!this.lengthKnown) {
                  this.cacheLength();
                }
                return this._chunkSize;
              },
            },
          });
          var properties = { isDevice: false, contents: lazyArray };
        } else {
          var properties = { isDevice: false, url: url };
        }
        var node = FS.createFile(parent, name, properties, canRead, canWrite);
        if (properties.contents) {
          node.contents = properties.contents;
        } else if (properties.url) {
          node.contents = null;
          node.url = properties.url;
        }
        Object.defineProperties(node, {
          usedBytes: {
            get: function () {
              return this.contents.length;
            },
          },
        });
        var stream_ops = {};
        var keys = Object.keys(node.stream_ops);
        keys.forEach((key) => {
          var fn = node.stream_ops[key];
          stream_ops[key] = function forceLoadLazyFile() {
            FS.forceLoadFile(node);
            return fn.apply(null, arguments);
          };
        });
        function writeChunks(stream, buffer, offset, length, position) {
          var contents = stream.node.contents;
          if (position >= contents.length) return 0;
          var size = Math.min(contents.length - position, length);
          assert(size >= 0);
          if (contents.slice) {
            for (var i = 0; i < size; i++) {
              buffer[offset + i] = contents[position + i];
            }
          } else {
            for (var i = 0; i < size; i++) {
              buffer[offset + i] = contents.get(position + i);
            }
          }
          return size;
        }
        stream_ops.read = (stream, buffer, offset, length, position) => {
          FS.forceLoadFile(node);
          return writeChunks(stream, buffer, offset, length, position);
        };
        stream_ops.mmap = (stream, length, position, prot, flags) => {
          FS.forceLoadFile(node);
          var ptr = mmapAlloc(length);
          if (!ptr) {
            throw new FS.ErrnoError(48);
          }
          writeChunks(stream, HEAP8, ptr, length, position);
          return { ptr: ptr, allocated: true };
        };
        node.stream_ops = stream_ops;
        return node;
      },
      absolutePath: () => {
        abort("FS.absolutePath has been removed; use PATH_FS.resolve instead");
      },
      createFolder: () => {
        abort("FS.createFolder has been removed; use FS.mkdir instead");
      },
      createLink: () => {
        abort("FS.createLink has been removed; use FS.symlink instead");
      },
      joinPath: () => {
        abort("FS.joinPath has been removed; use PATH.join instead");
      },
      mmapAlloc: () => {
        abort(
          "FS.mmapAlloc has been replaced by the top level function mmapAlloc",
        );
      },
      standardizePath: () => {
        abort(
          "FS.standardizePath has been removed; use PATH.normalize instead",
        );
      },
    };
    var SYSCALLS = {
      DEFAULT_POLLMASK: 5,
      calculateAt: function (dirfd, path, allowEmpty) {
        if (PATH.isAbs(path)) {
          return path;
        }
        var dir;
        if (dirfd === -100) {
          dir = FS.cwd();
        } else {
          var dirstream = SYSCALLS.getStreamFromFD(dirfd);
          dir = dirstream.path;
        }
        if (path.length == 0) {
          if (!allowEmpty) {
            throw new FS.ErrnoError(44);
          }
          return dir;
        }
        return PATH.join2(dir, path);
      },
      doStat: function (func, path, buf) {
        try {
          var stat = func(path);
        } catch (e) {
          if (
            e &&
            e.node &&
            PATH.normalize(path) !== PATH.normalize(FS.getPath(e.node))
          ) {
            return -54;
          }
          throw e;
        }
        HEAP32[buf >> 2] = stat.dev;
        checkInt32(stat.dev);
        HEAP32[(buf + 4) >> 2] = stat.mode;
        checkInt32(stat.mode);
        HEAPU32[(buf + 8) >> 2] = stat.nlink;
        checkInt32(stat.nlink);
        HEAP32[(buf + 12) >> 2] = stat.uid;
        checkInt32(stat.uid);
        HEAP32[(buf + 16) >> 2] = stat.gid;
        checkInt32(stat.gid);
        HEAP32[(buf + 20) >> 2] = stat.rdev;
        checkInt32(stat.rdev);
        (tempI64 = [
          stat.size >>> 0,
          ((tempDouble = stat.size),
          +Math.abs(tempDouble) >= 1
            ? tempDouble > 0
              ? +Math.floor(tempDouble / 4294967296) >>> 0
              : ~~+Math.ceil(
                  (tempDouble - +(~~tempDouble >>> 0)) / 4294967296,
                ) >>> 0
            : 0),
        ]),
          (HEAP32[(buf + 24) >> 2] = tempI64[0]),
          (HEAP32[(buf + 28) >> 2] = tempI64[1]);
        checkInt64(stat.size);
        HEAP32[(buf + 32) >> 2] = 4096;
        checkInt32(4096);
        HEAP32[(buf + 36) >> 2] = stat.blocks;
        checkInt32(stat.blocks);
        var atime = stat.atime.getTime();
        var mtime = stat.mtime.getTime();
        var ctime = stat.ctime.getTime();
        (tempI64 = [
          Math.floor(atime / 1e3) >>> 0,
          ((tempDouble = Math.floor(atime / 1e3)),
          +Math.abs(tempDouble) >= 1
            ? tempDouble > 0
              ? +Math.floor(tempDouble / 4294967296) >>> 0
              : ~~+Math.ceil(
                  (tempDouble - +(~~tempDouble >>> 0)) / 4294967296,
                ) >>> 0
            : 0),
        ]),
          (HEAP32[(buf + 40) >> 2] = tempI64[0]),
          (HEAP32[(buf + 44) >> 2] = tempI64[1]);
        checkInt64(Math.floor(atime / 1e3));
        HEAPU32[(buf + 48) >> 2] = (atime % 1e3) * 1e3;
        checkInt32((atime % 1e3) * 1e3);
        (tempI64 = [
          Math.floor(mtime / 1e3) >>> 0,
          ((tempDouble = Math.floor(mtime / 1e3)),
          +Math.abs(tempDouble) >= 1
            ? tempDouble > 0
              ? +Math.floor(tempDouble / 4294967296) >>> 0
              : ~~+Math.ceil(
                  (tempDouble - +(~~tempDouble >>> 0)) / 4294967296,
                ) >>> 0
            : 0),
        ]),
          (HEAP32[(buf + 56) >> 2] = tempI64[0]),
          (HEAP32[(buf + 60) >> 2] = tempI64[1]);
        checkInt64(Math.floor(mtime / 1e3));
        HEAPU32[(buf + 64) >> 2] = (mtime % 1e3) * 1e3;
        checkInt32((mtime % 1e3) * 1e3);
        (tempI64 = [
          Math.floor(ctime / 1e3) >>> 0,
          ((tempDouble = Math.floor(ctime / 1e3)),
          +Math.abs(tempDouble) >= 1
            ? tempDouble > 0
              ? +Math.floor(tempDouble / 4294967296) >>> 0
              : ~~+Math.ceil(
                  (tempDouble - +(~~tempDouble >>> 0)) / 4294967296,
                ) >>> 0
            : 0),
        ]),
          (HEAP32[(buf + 72) >> 2] = tempI64[0]),
          (HEAP32[(buf + 76) >> 2] = tempI64[1]);
        checkInt64(Math.floor(ctime / 1e3));
        HEAPU32[(buf + 80) >> 2] = (ctime % 1e3) * 1e3;
        checkInt32((ctime % 1e3) * 1e3);
        (tempI64 = [
          stat.ino >>> 0,
          ((tempDouble = stat.ino),
          +Math.abs(tempDouble) >= 1
            ? tempDouble > 0
              ? +Math.floor(tempDouble / 4294967296) >>> 0
              : ~~+Math.ceil(
                  (tempDouble - +(~~tempDouble >>> 0)) / 4294967296,
                ) >>> 0
            : 0),
        ]),
          (HEAP32[(buf + 88) >> 2] = tempI64[0]),
          (HEAP32[(buf + 92) >> 2] = tempI64[1]);
        checkInt64(stat.ino);
        return 0;
      },
      doMsync: function (addr, stream, len, flags, offset) {
        if (!FS.isFile(stream.node.mode)) {
          throw new FS.ErrnoError(43);
        }
        if (flags & 2) {
          return 0;
        }
        var buffer = HEAPU8.slice(addr, addr + len);
        FS.msync(stream, buffer, offset, len, flags);
      },
      varargs: undefined,
      get() {
        assert(SYSCALLS.varargs != undefined);
        SYSCALLS.varargs += 4;
        var ret = HEAP32[(SYSCALLS.varargs - 4) >> 2];
        return ret;
      },
      getStr(ptr) {
        var ret = UTF8ToString(ptr);
        return ret;
      },
      getStreamFromFD: function (fd) {
        var stream = FS.getStreamChecked(fd);
        return stream;
      },
    };
    function ___syscall_fcntl64(fd, cmd, varargs) {
      SYSCALLS.varargs = varargs;
      try {
        var stream = SYSCALLS.getStreamFromFD(fd);
        switch (cmd) {
          case 0: {
            var arg = SYSCALLS.get();
            if (arg < 0) {
              return -28;
            }
            var newStream;
            newStream = FS.createStream(stream, arg);
            return newStream.fd;
          }
          case 1:
          case 2:
            return 0;
          case 3:
            return stream.flags;
          case 4: {
            var arg = SYSCALLS.get();
            stream.flags |= arg;
            return 0;
          }
          case 5: {
            var arg = SYSCALLS.get();
            var offset = 0;
            HEAP16[(arg + offset) >> 1] = 2;
            checkInt16(2);
            return 0;
          }
          case 6:
          case 7:
            return 0;
          case 16:
          case 8:
            return -28;
          case 9:
            setErrNo(28);
            return -1;
          default: {
            return -28;
          }
        }
      } catch (e) {
        if (typeof FS == "undefined" || !(e.name === "ErrnoError")) throw e;
        return -e.errno;
      }
    }
    function ___syscall_ioctl(fd, op, varargs) {
      SYSCALLS.varargs = varargs;
      try {
        var stream = SYSCALLS.getStreamFromFD(fd);
        switch (op) {
          case 21509: {
            if (!stream.tty) return -59;
            return 0;
          }
          case 21505: {
            if (!stream.tty) return -59;
            if (stream.tty.ops.ioctl_tcgets) {
              var termios = stream.tty.ops.ioctl_tcgets(stream);
              var argp = SYSCALLS.get();
              HEAP32[argp >> 2] = termios.c_iflag || 0;
              checkInt32(termios.c_iflag || 0);
              HEAP32[(argp + 4) >> 2] = termios.c_oflag || 0;
              checkInt32(termios.c_oflag || 0);
              HEAP32[(argp + 8) >> 2] = termios.c_cflag || 0;
              checkInt32(termios.c_cflag || 0);
              HEAP32[(argp + 12) >> 2] = termios.c_lflag || 0;
              checkInt32(termios.c_lflag || 0);
              for (var i = 0; i < 32; i++) {
                HEAP8[(argp + i + 17) >> 0] = termios.c_cc[i] || 0;
                checkInt8(termios.c_cc[i] || 0);
              }
              return 0;
            }
            return 0;
          }
          case 21510:
          case 21511:
          case 21512: {
            if (!stream.tty) return -59;
            return 0;
          }
          case 21506:
          case 21507:
          case 21508: {
            if (!stream.tty) return -59;
            if (stream.tty.ops.ioctl_tcsets) {
              var argp = SYSCALLS.get();
              var c_iflag = HEAP32[argp >> 2];
              var c_oflag = HEAP32[(argp + 4) >> 2];
              var c_cflag = HEAP32[(argp + 8) >> 2];
              var c_lflag = HEAP32[(argp + 12) >> 2];
              var c_cc = [];
              for (var i = 0; i < 32; i++) {
                c_cc.push(HEAP8[(argp + i + 17) >> 0]);
              }
              return stream.tty.ops.ioctl_tcsets(stream.tty, op, {
                c_iflag: c_iflag,
                c_oflag: c_oflag,
                c_cflag: c_cflag,
                c_lflag: c_lflag,
                c_cc: c_cc,
              });
            }
            return 0;
          }
          case 21519: {
            if (!stream.tty) return -59;
            var argp = SYSCALLS.get();
            HEAP32[argp >> 2] = 0;
            checkInt32(0);
            return 0;
          }
          case 21520: {
            if (!stream.tty) return -59;
            return -28;
          }
          case 21531: {
            var argp = SYSCALLS.get();
            return FS.ioctl(stream, op, argp);
          }
          case 21523: {
            if (!stream.tty) return -59;
            if (stream.tty.ops.ioctl_tiocgwinsz) {
              var winsize = stream.tty.ops.ioctl_tiocgwinsz(stream.tty);
              var argp = SYSCALLS.get();
              HEAP16[argp >> 1] = winsize[0];
              checkInt16(winsize[0]);
              HEAP16[(argp + 2) >> 1] = winsize[1];
              checkInt16(winsize[1]);
            }
            return 0;
          }
          case 21524: {
            if (!stream.tty) return -59;
            return 0;
          }
          case 21515: {
            if (!stream.tty) return -59;
            return 0;
          }
          default:
            return -28;
        }
      } catch (e) {
        if (typeof FS == "undefined" || !(e.name === "ErrnoError")) throw e;
        return -e.errno;
      }
    }
    function ___syscall_openat(dirfd, path, flags, varargs) {
      SYSCALLS.varargs = varargs;
      try {
        path = SYSCALLS.getStr(path);
        path = SYSCALLS.calculateAt(dirfd, path);
        var mode = varargs ? SYSCALLS.get() : 0;
        return FS.open(path, flags, mode).fd;
      } catch (e) {
        if (typeof FS == "undefined" || !(e.name === "ErrnoError")) throw e;
        return -e.errno;
      }
    }
    function ___throw_exception_with_stack_trace(ex) {
      var e = new WebAssembly.Exception(getCppExceptionTag(), [ex], {
        traceStack: true,
      });
      e.message = getExceptionMessage(e);
      if (e.stack) {
        var arr = e.stack.split("\n");
        arr.splice(1, 1);
        e.stack = arr.join("\n");
      }
      throw e;
    }
    var structRegistrations = {};
    function runDestructors(destructors) {
      while (destructors.length) {
        var ptr = destructors.pop();
        var del = destructors.pop();
        del(ptr);
      }
    }
    function simpleReadValueFromPointer(pointer) {
      return this["fromWireType"](HEAP32[pointer >> 2]);
    }
    var awaitingDependencies = {};
    var registeredTypes = {};
    var typeDependencies = {};
    var InternalError = undefined;
    function throwInternalError(message) {
      throw new InternalError(message);
    }
    function whenDependentTypesAreResolved(
      myTypes,
      dependentTypes,
      getTypeConverters,
    ) {
      myTypes.forEach(function (type) {
        typeDependencies[type] = dependentTypes;
      });
      function onComplete(typeConverters) {
        var myTypeConverters = getTypeConverters(typeConverters);
        if (myTypeConverters.length !== myTypes.length) {
          throwInternalError("Mismatched type converter count");
        }
        for (var i = 0; i < myTypes.length; ++i) {
          registerType(myTypes[i], myTypeConverters[i]);
        }
      }
      var typeConverters = new Array(dependentTypes.length);
      var unregisteredTypes = [];
      var registered = 0;
      dependentTypes.forEach((dt, i) => {
        if (registeredTypes.hasOwnProperty(dt)) {
          typeConverters[i] = registeredTypes[dt];
        } else {
          unregisteredTypes.push(dt);
          if (!awaitingDependencies.hasOwnProperty(dt)) {
            awaitingDependencies[dt] = [];
          }
          awaitingDependencies[dt].push(() => {
            typeConverters[i] = registeredTypes[dt];
            ++registered;
            if (registered === unregisteredTypes.length) {
              onComplete(typeConverters);
            }
          });
        }
      });
      if (0 === unregisteredTypes.length) {
        onComplete(typeConverters);
      }
    }
    var __embind_finalize_value_object = function (structType) {
      var reg = structRegistrations[structType];
      delete structRegistrations[structType];
      var rawConstructor = reg.rawConstructor;
      var rawDestructor = reg.rawDestructor;
      var fieldRecords = reg.fields;
      var fieldTypes = fieldRecords
        .map((field) => field.getterReturnType)
        .concat(fieldRecords.map((field) => field.setterArgumentType));
      whenDependentTypesAreResolved([structType], fieldTypes, (fieldTypes) => {
        var fields = {};
        fieldRecords.forEach((field, i) => {
          var fieldName = field.fieldName;
          var getterReturnType = fieldTypes[i];
          var getter = field.getter;
          var getterContext = field.getterContext;
          var setterArgumentType = fieldTypes[i + fieldRecords.length];
          var setter = field.setter;
          var setterContext = field.setterContext;
          fields[fieldName] = {
            read: (ptr) =>
              getterReturnType["fromWireType"](getter(getterContext, ptr)),
            write: (ptr, o) => {
              var destructors = [];
              setter(
                setterContext,
                ptr,
                setterArgumentType["toWireType"](destructors, o),
              );
              runDestructors(destructors);
            },
          };
        });
        return [
          {
            name: reg.name,
            fromWireType: function (ptr) {
              var rv = {};
              for (var i in fields) {
                rv[i] = fields[i].read(ptr);
              }
              rawDestructor(ptr);
              return rv;
            },
            toWireType: function (destructors, o) {
              for (var fieldName in fields) {
                if (!(fieldName in o)) {
                  throw new TypeError(`Missing field: "${fieldName}"`);
                }
              }
              var ptr = rawConstructor();
              for (fieldName in fields) {
                fields[fieldName].write(ptr, o[fieldName]);
              }
              if (destructors !== null) {
                destructors.push(rawDestructor, ptr);
              }
              return ptr;
            },
            argPackAdvance: 8,
            readValueFromPointer: simpleReadValueFromPointer,
            destructorFunction: rawDestructor,
          },
        ];
      });
    };
    function __embind_register_bigint(
      primitiveType,
      name,
      size,
      minRange,
      maxRange,
    ) {}
    function getShiftFromSize(size) {
      switch (size) {
        case 1:
          return 0;
        case 2:
          return 1;
        case 4:
          return 2;
        case 8:
          return 3;
        default:
          throw new TypeError(`Unknown type size: ${size}`);
      }
    }
    function embind_init_charCodes() {
      var codes = new Array(256);
      for (var i = 0; i < 256; ++i) {
        codes[i] = String.fromCharCode(i);
      }
      embind_charCodes = codes;
    }
    var embind_charCodes = undefined;
    function readLatin1String(ptr) {
      var ret = "";
      var c = ptr;
      while (HEAPU8[c]) {
        ret += embind_charCodes[HEAPU8[c++]];
      }
      return ret;
    }
    var BindingError = undefined;
    function throwBindingError(message) {
      throw new BindingError(message);
    }
    function sharedRegisterType(rawType, registeredInstance, options = {}) {
      var name = registeredInstance.name;
      if (!rawType) {
        throwBindingError(
          `type "${name}" must have a positive integer typeid pointer`,
        );
      }
      if (registeredTypes.hasOwnProperty(rawType)) {
        if (options.ignoreDuplicateRegistrations) {
          return;
        } else {
          throwBindingError(`Cannot register type '${name}' twice`);
        }
      }
      registeredTypes[rawType] = registeredInstance;
      delete typeDependencies[rawType];
      if (awaitingDependencies.hasOwnProperty(rawType)) {
        var callbacks = awaitingDependencies[rawType];
        delete awaitingDependencies[rawType];
        callbacks.forEach((cb) => cb());
      }
    }
    function registerType(rawType, registeredInstance, options = {}) {
      if (!("argPackAdvance" in registeredInstance)) {
        throw new TypeError(
          "registerType registeredInstance requires argPackAdvance",
        );
      }
      return sharedRegisterType(rawType, registeredInstance, options);
    }
    function __embind_register_bool(
      rawType,
      name,
      size,
      trueValue,
      falseValue,
    ) {
      var shift = getShiftFromSize(size);
      name = readLatin1String(name);
      registerType(rawType, {
        name: name,
        fromWireType: function (wt) {
          return !!wt;
        },
        toWireType: function (destructors, o) {
          return o ? trueValue : falseValue;
        },
        argPackAdvance: 8,
        readValueFromPointer: function (pointer) {
          var heap;
          if (size === 1) {
            heap = HEAP8;
          } else if (size === 2) {
            heap = HEAP16;
          } else if (size === 4) {
            heap = HEAP32;
          } else {
            throw new TypeError("Unknown boolean type size: " + name);
          }
          return this["fromWireType"](heap[pointer >> shift]);
        },
        destructorFunction: null,
      });
    }
    function ClassHandle_isAliasOf(other) {
      if (!(this instanceof ClassHandle)) {
        return false;
      }
      if (!(other instanceof ClassHandle)) {
        return false;
      }
      var leftClass = this.$$.ptrType.registeredClass;
      var left = this.$$.ptr;
      var rightClass = other.$$.ptrType.registeredClass;
      var right = other.$$.ptr;
      while (leftClass.baseClass) {
        left = leftClass.upcast(left);
        leftClass = leftClass.baseClass;
      }
      while (rightClass.baseClass) {
        right = rightClass.upcast(right);
        rightClass = rightClass.baseClass;
      }
      return leftClass === rightClass && left === right;
    }
    function shallowCopyInternalPointer(o) {
      return {
        count: o.count,
        deleteScheduled: o.deleteScheduled,
        preservePointerOnDelete: o.preservePointerOnDelete,
        ptr: o.ptr,
        ptrType: o.ptrType,
        smartPtr: o.smartPtr,
        smartPtrType: o.smartPtrType,
      };
    }
    function throwInstanceAlreadyDeleted(obj) {
      function getInstanceTypeName(handle) {
        return handle.$$.ptrType.registeredClass.name;
      }
      throwBindingError(getInstanceTypeName(obj) + " instance already deleted");
    }
    var finalizationRegistry = false;
    function detachFinalizer(handle) {}
    function runDestructor($$) {
      if ($$.smartPtr) {
        $$.smartPtrType.rawDestructor($$.smartPtr);
      } else {
        $$.ptrType.registeredClass.rawDestructor($$.ptr);
      }
    }
    function releaseClassHandle($$) {
      $$.count.value -= 1;
      var toDelete = 0 === $$.count.value;
      if (toDelete) {
        runDestructor($$);
      }
    }
    function downcastPointer(ptr, ptrClass, desiredClass) {
      if (ptrClass === desiredClass) {
        return ptr;
      }
      if (undefined === desiredClass.baseClass) {
        return null;
      }
      var rv = downcastPointer(ptr, ptrClass, desiredClass.baseClass);
      if (rv === null) {
        return null;
      }
      return desiredClass.downcast(rv);
    }
    var registeredPointers = {};
    function getInheritedInstanceCount() {
      return Object.keys(registeredInstances).length;
    }
    function getLiveInheritedInstances() {
      var rv = [];
      for (var k in registeredInstances) {
        if (registeredInstances.hasOwnProperty(k)) {
          rv.push(registeredInstances[k]);
        }
      }
      return rv;
    }
    var deletionQueue = [];
    function flushPendingDeletes() {
      while (deletionQueue.length) {
        var obj = deletionQueue.pop();
        obj.$$.deleteScheduled = false;
        obj["delete"]();
      }
    }
    var delayFunction = undefined;
    function setDelayFunction(fn) {
      delayFunction = fn;
      if (deletionQueue.length && delayFunction) {
        delayFunction(flushPendingDeletes);
      }
    }
    function init_embind() {
      Module["getInheritedInstanceCount"] = getInheritedInstanceCount;
      Module["getLiveInheritedInstances"] = getLiveInheritedInstances;
      Module["flushPendingDeletes"] = flushPendingDeletes;
      Module["setDelayFunction"] = setDelayFunction;
    }
    var registeredInstances = {};
    function getBasestPointer(class_, ptr) {
      if (ptr === undefined) {
        throwBindingError("ptr should not be undefined");
      }
      while (class_.baseClass) {
        ptr = class_.upcast(ptr);
        class_ = class_.baseClass;
      }
      return ptr;
    }
    function getInheritedInstance(class_, ptr) {
      ptr = getBasestPointer(class_, ptr);
      return registeredInstances[ptr];
    }
    function makeClassHandle(prototype, record) {
      if (!record.ptrType || !record.ptr) {
        throwInternalError("makeClassHandle requires ptr and ptrType");
      }
      var hasSmartPtrType = !!record.smartPtrType;
      var hasSmartPtr = !!record.smartPtr;
      if (hasSmartPtrType !== hasSmartPtr) {
        throwInternalError("Both smartPtrType and smartPtr must be specified");
      }
      record.count = { value: 1 };
      return attachFinalizer(
        Object.create(prototype, { $$: { value: record } }),
      );
    }
    function RegisteredPointer_fromWireType(ptr) {
      var rawPointer = this.getPointee(ptr);
      if (!rawPointer) {
        this.destructor(ptr);
        return null;
      }
      var registeredInstance = getInheritedInstance(
        this.registeredClass,
        rawPointer,
      );
      if (undefined !== registeredInstance) {
        if (0 === registeredInstance.$$.count.value) {
          registeredInstance.$$.ptr = rawPointer;
          registeredInstance.$$.smartPtr = ptr;
          return registeredInstance["clone"]();
        } else {
          var rv = registeredInstance["clone"]();
          this.destructor(ptr);
          return rv;
        }
      }
      function makeDefaultHandle() {
        if (this.isSmartPointer) {
          return makeClassHandle(this.registeredClass.instancePrototype, {
            ptrType: this.pointeeType,
            ptr: rawPointer,
            smartPtrType: this,
            smartPtr: ptr,
          });
        } else {
          return makeClassHandle(this.registeredClass.instancePrototype, {
            ptrType: this,
            ptr: ptr,
          });
        }
      }
      var actualType = this.registeredClass.getActualType(rawPointer);
      var registeredPointerRecord = registeredPointers[actualType];
      if (!registeredPointerRecord) {
        return makeDefaultHandle.call(this);
      }
      var toType;
      if (this.isConst) {
        toType = registeredPointerRecord.constPointerType;
      } else {
        toType = registeredPointerRecord.pointerType;
      }
      var dp = downcastPointer(
        rawPointer,
        this.registeredClass,
        toType.registeredClass,
      );
      if (dp === null) {
        return makeDefaultHandle.call(this);
      }
      if (this.isSmartPointer) {
        return makeClassHandle(toType.registeredClass.instancePrototype, {
          ptrType: toType,
          ptr: dp,
          smartPtrType: this,
          smartPtr: ptr,
        });
      } else {
        return makeClassHandle(toType.registeredClass.instancePrototype, {
          ptrType: toType,
          ptr: dp,
        });
      }
    }
    var attachFinalizer = function (handle) {
      if ("undefined" === typeof FinalizationRegistry) {
        attachFinalizer = (handle) => handle;
        return handle;
      }
      finalizationRegistry = new FinalizationRegistry((info) => {
        console.warn(info.leakWarning.stack.replace(/^Error: /, ""));
        releaseClassHandle(info.$$);
      });
      attachFinalizer = (handle) => {
        var $$ = handle.$$;
        var hasSmartPtr = !!$$.smartPtr;
        if (hasSmartPtr) {
          var info = { $$: $$ };
          var cls = $$.ptrType.registeredClass;
          info.leakWarning = new Error(
            `Embind found a leaked C++ instance ${cls.name} <${ptrToString(
              $$.ptr,
            )}>.\n` +
              "We'll free it automatically in this case, but this functionality is not reliable across various environments.\n" +
              "Make sure to invoke .delete() manually once you're done with the instance instead.\n" +
              "Originally allocated",
          );
          if ("captureStackTrace" in Error) {
            Error.captureStackTrace(
              info.leakWarning,
              RegisteredPointer_fromWireType,
            );
          }
          finalizationRegistry.register(handle, info, handle);
        }
        return handle;
      };
      detachFinalizer = (handle) => finalizationRegistry.unregister(handle);
      return attachFinalizer(handle);
    };
    function ClassHandle_clone() {
      if (!this.$$.ptr) {
        throwInstanceAlreadyDeleted(this);
      }
      if (this.$$.preservePointerOnDelete) {
        this.$$.count.value += 1;
        return this;
      } else {
        var clone = attachFinalizer(
          Object.create(Object.getPrototypeOf(this), {
            $$: { value: shallowCopyInternalPointer(this.$$) },
          }),
        );
        clone.$$.count.value += 1;
        clone.$$.deleteScheduled = false;
        return clone;
      }
    }
    function ClassHandle_delete() {
      if (!this.$$.ptr) {
        throwInstanceAlreadyDeleted(this);
      }
      if (this.$$.deleteScheduled && !this.$$.preservePointerOnDelete) {
        throwBindingError("Object already scheduled for deletion");
      }
      detachFinalizer(this);
      releaseClassHandle(this.$$);
      if (!this.$$.preservePointerOnDelete) {
        this.$$.smartPtr = undefined;
        this.$$.ptr = undefined;
      }
    }
    function ClassHandle_isDeleted() {
      return !this.$$.ptr;
    }
    function ClassHandle_deleteLater() {
      if (!this.$$.ptr) {
        throwInstanceAlreadyDeleted(this);
      }
      if (this.$$.deleteScheduled && !this.$$.preservePointerOnDelete) {
        throwBindingError("Object already scheduled for deletion");
      }
      deletionQueue.push(this);
      if (deletionQueue.length === 1 && delayFunction) {
        delayFunction(flushPendingDeletes);
      }
      this.$$.deleteScheduled = true;
      return this;
    }
    function init_ClassHandle() {
      ClassHandle.prototype["isAliasOf"] = ClassHandle_isAliasOf;
      ClassHandle.prototype["clone"] = ClassHandle_clone;
      ClassHandle.prototype["delete"] = ClassHandle_delete;
      ClassHandle.prototype["isDeleted"] = ClassHandle_isDeleted;
      ClassHandle.prototype["deleteLater"] = ClassHandle_deleteLater;
    }
    function ClassHandle() {}
    var char_0 = 48;
    var char_9 = 57;
    function makeLegalFunctionName(name) {
      if (undefined === name) {
        return "_unknown";
      }
      name = name.replace(/[^a-zA-Z0-9_]/g, "$");
      var f = name.charCodeAt(0);
      if (f >= char_0 && f <= char_9) {
        return `_${name}`;
      }
      return name;
    }
    function createNamedFunction(name, body) {
      name = makeLegalFunctionName(name);
      return {
        [name]: function () {
          return body.apply(this, arguments);
        },
      }[name];
    }
    function ensureOverloadTable(proto, methodName, humanName) {
      if (undefined === proto[methodName].overloadTable) {
        var prevFunc = proto[methodName];
        proto[methodName] = function () {
          if (
            !proto[methodName].overloadTable.hasOwnProperty(arguments.length)
          ) {
            throwBindingError(
              `Function '${humanName}' called with an invalid number of arguments (${arguments.length}) - expects one of (${proto[methodName].overloadTable})!`,
            );
          }
          return proto[methodName].overloadTable[arguments.length].apply(
            this,
            arguments,
          );
        };
        proto[methodName].overloadTable = [];
        proto[methodName].overloadTable[prevFunc.argCount] = prevFunc;
      }
    }
    function exposePublicSymbol(name, value, numArguments) {
      if (Module.hasOwnProperty(name)) {
        if (
          undefined === numArguments ||
          (undefined !== Module[name].overloadTable &&
            undefined !== Module[name].overloadTable[numArguments])
        ) {
          throwBindingError(`Cannot register public name '${name}' twice`);
        }
        ensureOverloadTable(Module, name, name);
        if (Module.hasOwnProperty(numArguments)) {
          throwBindingError(
            `Cannot register multiple overloads of a function with the same number of arguments (${numArguments})!`,
          );
        }
        Module[name].overloadTable[numArguments] = value;
      } else {
        Module[name] = value;
        if (undefined !== numArguments) {
          Module[name].numArguments = numArguments;
        }
      }
    }
    function RegisteredClass(
      name,
      constructor,
      instancePrototype,
      rawDestructor,
      baseClass,
      getActualType,
      upcast,
      downcast,
    ) {
      this.name = name;
      this.constructor = constructor;
      this.instancePrototype = instancePrototype;
      this.rawDestructor = rawDestructor;
      this.baseClass = baseClass;
      this.getActualType = getActualType;
      this.upcast = upcast;
      this.downcast = downcast;
      this.pureVirtualFunctions = [];
    }
    function upcastPointer(ptr, ptrClass, desiredClass) {
      while (ptrClass !== desiredClass) {
        if (!ptrClass.upcast) {
          throwBindingError(
            `Expected null or instance of ${desiredClass.name}, got an instance of ${ptrClass.name}`,
          );
        }
        ptr = ptrClass.upcast(ptr);
        ptrClass = ptrClass.baseClass;
      }
      return ptr;
    }
    function constNoSmartPtrRawPointerToWireType(destructors, handle) {
      if (handle === null) {
        if (this.isReference) {
          throwBindingError(`null is not a valid ${this.name}`);
        }
        return 0;
      }
      if (!handle.$$) {
        throwBindingError(
          `Cannot pass "${embindRepr(handle)}" as a ${this.name}`,
        );
      }
      if (!handle.$$.ptr) {
        throwBindingError(
          `Cannot pass deleted object as a pointer of type ${this.name}`,
        );
      }
      var handleClass = handle.$$.ptrType.registeredClass;
      var ptr = upcastPointer(handle.$$.ptr, handleClass, this.registeredClass);
      return ptr;
    }
    function genericPointerToWireType(destructors, handle) {
      var ptr;
      if (handle === null) {
        if (this.isReference) {
          throwBindingError(`null is not a valid ${this.name}`);
        }
        if (this.isSmartPointer) {
          ptr = this.rawConstructor();
          if (destructors !== null) {
            destructors.push(this.rawDestructor, ptr);
          }
          return ptr;
        } else {
          return 0;
        }
      }
      if (!handle.$$) {
        throwBindingError(
          `Cannot pass "${embindRepr(handle)}" as a ${this.name}`,
        );
      }
      if (!handle.$$.ptr) {
        throwBindingError(
          `Cannot pass deleted object as a pointer of type ${this.name}`,
        );
      }
      if (!this.isConst && handle.$$.ptrType.isConst) {
        throwBindingError(
          `Cannot convert argument of type ${
            handle.$$.smartPtrType
              ? handle.$$.smartPtrType.name
              : handle.$$.ptrType.name
          } to parameter type ${this.name}`,
        );
      }
      var handleClass = handle.$$.ptrType.registeredClass;
      ptr = upcastPointer(handle.$$.ptr, handleClass, this.registeredClass);
      if (this.isSmartPointer) {
        if (undefined === handle.$$.smartPtr) {
          throwBindingError("Passing raw pointer to smart pointer is illegal");
        }
        switch (this.sharingPolicy) {
          case 0:
            if (handle.$$.smartPtrType === this) {
              ptr = handle.$$.smartPtr;
            } else {
              throwBindingError(
                `Cannot convert argument of type ${
                  handle.$$.smartPtrType
                    ? handle.$$.smartPtrType.name
                    : handle.$$.ptrType.name
                } to parameter type ${this.name}`,
              );
            }
            break;
          case 1:
            ptr = handle.$$.smartPtr;
            break;
          case 2:
            if (handle.$$.smartPtrType === this) {
              ptr = handle.$$.smartPtr;
            } else {
              var clonedHandle = handle["clone"]();
              ptr = this.rawShare(
                ptr,
                Emval.toHandle(function () {
                  clonedHandle["delete"]();
                }),
              );
              if (destructors !== null) {
                destructors.push(this.rawDestructor, ptr);
              }
            }
            break;
          default:
            throwBindingError("Unsupporting sharing policy");
        }
      }
      return ptr;
    }
    function nonConstNoSmartPtrRawPointerToWireType(destructors, handle) {
      if (handle === null) {
        if (this.isReference) {
          throwBindingError(`null is not a valid ${this.name}`);
        }
        return 0;
      }
      if (!handle.$$) {
        throwBindingError(
          `Cannot pass "${embindRepr(handle)}" as a ${this.name}`,
        );
      }
      if (!handle.$$.ptr) {
        throwBindingError(
          `Cannot pass deleted object as a pointer of type ${this.name}`,
        );
      }
      if (handle.$$.ptrType.isConst) {
        throwBindingError(
          `Cannot convert argument of type ${handle.$$.ptrType.name} to parameter type ${this.name}`,
        );
      }
      var handleClass = handle.$$.ptrType.registeredClass;
      var ptr = upcastPointer(handle.$$.ptr, handleClass, this.registeredClass);
      return ptr;
    }
    function RegisteredPointer_getPointee(ptr) {
      if (this.rawGetPointee) {
        ptr = this.rawGetPointee(ptr);
      }
      return ptr;
    }
    function RegisteredPointer_destructor(ptr) {
      if (this.rawDestructor) {
        this.rawDestructor(ptr);
      }
    }
    function RegisteredPointer_deleteObject(handle) {
      if (handle !== null) {
        handle["delete"]();
      }
    }
    function init_RegisteredPointer() {
      RegisteredPointer.prototype.getPointee = RegisteredPointer_getPointee;
      RegisteredPointer.prototype.destructor = RegisteredPointer_destructor;
      RegisteredPointer.prototype["argPackAdvance"] = 8;
      RegisteredPointer.prototype["readValueFromPointer"] =
        simpleReadValueFromPointer;
      RegisteredPointer.prototype["deleteObject"] =
        RegisteredPointer_deleteObject;
      RegisteredPointer.prototype["fromWireType"] =
        RegisteredPointer_fromWireType;
    }
    function RegisteredPointer(
      name,
      registeredClass,
      isReference,
      isConst,
      isSmartPointer,
      pointeeType,
      sharingPolicy,
      rawGetPointee,
      rawConstructor,
      rawShare,
      rawDestructor,
    ) {
      this.name = name;
      this.registeredClass = registeredClass;
      this.isReference = isReference;
      this.isConst = isConst;
      this.isSmartPointer = isSmartPointer;
      this.pointeeType = pointeeType;
      this.sharingPolicy = sharingPolicy;
      this.rawGetPointee = rawGetPointee;
      this.rawConstructor = rawConstructor;
      this.rawShare = rawShare;
      this.rawDestructor = rawDestructor;
      if (!isSmartPointer && registeredClass.baseClass === undefined) {
        if (isConst) {
          this["toWireType"] = constNoSmartPtrRawPointerToWireType;
          this.destructorFunction = null;
        } else {
          this["toWireType"] = nonConstNoSmartPtrRawPointerToWireType;
          this.destructorFunction = null;
        }
      } else {
        this["toWireType"] = genericPointerToWireType;
      }
    }
    function replacePublicSymbol(name, value, numArguments) {
      if (!Module.hasOwnProperty(name)) {
        throwInternalError("Replacing nonexistant public symbol");
      }
      if (
        undefined !== Module[name].overloadTable &&
        undefined !== numArguments
      ) {
        Module[name].overloadTable[numArguments] = value;
      } else {
        Module[name] = value;
        Module[name].argCount = numArguments;
      }
    }
    var dynCallLegacy = (sig, ptr, args) => {
      assert(
        "dynCall_" + sig in Module,
        `bad function pointer type - dynCall function not found for sig '${sig}'`,
      );
      if (args && args.length) {
        assert(args.length === sig.substring(1).replace(/j/g, "--").length);
      } else {
        assert(sig.length == 1);
      }
      var f = Module["dynCall_" + sig];
      return args && args.length
        ? f.apply(null, [ptr].concat(args))
        : f.call(null, ptr);
    };
    var wasmTableMirror = [];
    var getWasmTableEntry = (funcPtr) => {
      var func = wasmTableMirror[funcPtr];
      if (!func) {
        if (funcPtr >= wasmTableMirror.length)
          wasmTableMirror.length = funcPtr + 1;
        wasmTableMirror[funcPtr] = func = wasmTable.get(funcPtr);
      }
      assert(
        wasmTable.get(funcPtr) == func,
        "JavaScript-side Wasm function table mirror is out of date!",
      );
      return func;
    };
    var dynCall = (sig, ptr, args) => {
      if (sig.includes("j")) {
        return dynCallLegacy(sig, ptr, args);
      }
      assert(getWasmTableEntry(ptr), `missing table entry in dynCall: ${ptr}`);
      var rtn = getWasmTableEntry(ptr).apply(null, args);
      return rtn;
    };
    var getDynCaller = (sig, ptr) => {
      assert(
        sig.includes("j") || sig.includes("p"),
        "getDynCaller should only be called with i64 sigs",
      );
      var argCache = [];
      return function () {
        argCache.length = 0;
        Object.assign(argCache, arguments);
        return dynCall(sig, ptr, argCache);
      };
    };
    function embind__requireFunction(signature, rawFunction) {
      signature = readLatin1String(signature);
      function makeDynCaller() {
        if (signature.includes("j")) {
          return getDynCaller(signature, rawFunction);
        }
        return getWasmTableEntry(rawFunction);
      }
      var fp = makeDynCaller();
      if (typeof fp != "function") {
        throwBindingError(
          `unknown function pointer with signature ${signature}: ${rawFunction}`,
        );
      }
      return fp;
    }
    function extendError(baseErrorType, errorName) {
      var errorClass = createNamedFunction(errorName, function (message) {
        this.name = errorName;
        this.message = message;
        var stack = new Error(message).stack;
        if (stack !== undefined) {
          this.stack =
            this.toString() + "\n" + stack.replace(/^Error(:[^\n]*)?\n/, "");
        }
      });
      errorClass.prototype = Object.create(baseErrorType.prototype);
      errorClass.prototype.constructor = errorClass;
      errorClass.prototype.toString = function () {
        if (this.message === undefined) {
          return this.name;
        } else {
          return `${this.name}: ${this.message}`;
        }
      };
      return errorClass;
    }
    var UnboundTypeError = undefined;
    function getTypeName(type) {
      var ptr = ___getTypeName(type);
      var rv = readLatin1String(ptr);
      _free(ptr);
      return rv;
    }
    function throwUnboundTypeError(message, types) {
      var unboundTypes = [];
      var seen = {};
      function visit(type) {
        if (seen[type]) {
          return;
        }
        if (registeredTypes[type]) {
          return;
        }
        if (typeDependencies[type]) {
          typeDependencies[type].forEach(visit);
          return;
        }
        unboundTypes.push(type);
        seen[type] = true;
      }
      types.forEach(visit);
      throw new UnboundTypeError(
        `${message}: ` + unboundTypes.map(getTypeName).join([", "]),
      );
    }
    function __embind_register_class(
      rawType,
      rawPointerType,
      rawConstPointerType,
      baseClassRawType,
      getActualTypeSignature,
      getActualType,
      upcastSignature,
      upcast,
      downcastSignature,
      downcast,
      name,
      destructorSignature,
      rawDestructor,
    ) {
      name = readLatin1String(name);
      getActualType = embind__requireFunction(
        getActualTypeSignature,
        getActualType,
      );
      if (upcast) {
        upcast = embind__requireFunction(upcastSignature, upcast);
      }
      if (downcast) {
        downcast = embind__requireFunction(downcastSignature, downcast);
      }
      rawDestructor = embind__requireFunction(
        destructorSignature,
        rawDestructor,
      );
      var legalFunctionName = makeLegalFunctionName(name);
      exposePublicSymbol(legalFunctionName, function () {
        throwUnboundTypeError(`Cannot construct ${name} due to unbound types`, [
          baseClassRawType,
        ]);
      });
      whenDependentTypesAreResolved(
        [rawType, rawPointerType, rawConstPointerType],
        baseClassRawType ? [baseClassRawType] : [],
        function (base) {
          base = base[0];
          var baseClass;
          var basePrototype;
          if (baseClassRawType) {
            baseClass = base.registeredClass;
            basePrototype = baseClass.instancePrototype;
          } else {
            basePrototype = ClassHandle.prototype;
          }
          var constructor = createNamedFunction(legalFunctionName, function () {
            if (Object.getPrototypeOf(this) !== instancePrototype) {
              throw new BindingError("Use 'new' to construct " + name);
            }
            if (undefined === registeredClass.constructor_body) {
              throw new BindingError(name + " has no accessible constructor");
            }
            var body = registeredClass.constructor_body[arguments.length];
            if (undefined === body) {
              throw new BindingError(
                `Tried to invoke ctor of ${name} with invalid number of parameters (${
                  arguments.length
                }) - expected (${Object.keys(
                  registeredClass.constructor_body,
                ).toString()}) parameters instead!`,
              );
            }
            return body.apply(this, arguments);
          });
          var instancePrototype = Object.create(basePrototype, {
            constructor: { value: constructor },
          });
          constructor.prototype = instancePrototype;
          var registeredClass = new RegisteredClass(
            name,
            constructor,
            instancePrototype,
            rawDestructor,
            baseClass,
            getActualType,
            upcast,
            downcast,
          );
          if (registeredClass.baseClass) {
            if (registeredClass.baseClass.__derivedClasses === undefined) {
              registeredClass.baseClass.__derivedClasses = [];
            }
            registeredClass.baseClass.__derivedClasses.push(registeredClass);
          }
          var referenceConverter = new RegisteredPointer(
            name,
            registeredClass,
            true,
            false,
            false,
          );
          var pointerConverter = new RegisteredPointer(
            name + "*",
            registeredClass,
            false,
            false,
            false,
          );
          var constPointerConverter = new RegisteredPointer(
            name + " const*",
            registeredClass,
            false,
            true,
            false,
          );
          registeredPointers[rawType] = {
            pointerType: pointerConverter,
            constPointerType: constPointerConverter,
          };
          replacePublicSymbol(legalFunctionName, constructor);
          return [referenceConverter, pointerConverter, constPointerConverter];
        },
      );
    }
    function heap32VectorToArray(count, firstElement) {
      var array = [];
      for (var i = 0; i < count; i++) {
        array.push(HEAPU32[(firstElement + i * 4) >> 2]);
      }
      return array;
    }
    function newFunc(constructor, argumentList) {
      if (!(constructor instanceof Function)) {
        throw new TypeError(
          `new_ called with constructor type ${typeof constructor} which is not a function`,
        );
      }
      var dummy = createNamedFunction(
        constructor.name || "unknownFunctionName",
        function () {},
      );
      dummy.prototype = constructor.prototype;
      var obj = new dummy();
      var r = constructor.apply(obj, argumentList);
      return r instanceof Object ? r : obj;
    }
    function craftInvokerFunction(
      humanName,
      argTypes,
      classType,
      cppInvokerFunc,
      cppTargetFunc,
      isAsync,
    ) {
      var argCount = argTypes.length;
      if (argCount < 2) {
        throwBindingError(
          "argTypes array size mismatch! Must at least get return value and 'this' types!",
        );
      }
      assert(!isAsync, "Async bindings are only supported with JSPI.");
      var isClassMethodFunc = argTypes[1] !== null && classType !== null;
      var needsDestructorStack = false;
      for (var i = 1; i < argTypes.length; ++i) {
        if (
          argTypes[i] !== null &&
          argTypes[i].destructorFunction === undefined
        ) {
          needsDestructorStack = true;
          break;
        }
      }
      var returns = argTypes[0].name !== "void";
      var argsList = "";
      var argsListWired = "";
      for (var i = 0; i < argCount - 2; ++i) {
        argsList += (i !== 0 ? ", " : "") + "arg" + i;
        argsListWired += (i !== 0 ? ", " : "") + "arg" + i + "Wired";
      }
      var invokerFnBody = `\n        return function ${makeLegalFunctionName(
        humanName,
      )}(${argsList}) {\n        if (arguments.length !== ${
        argCount - 2
      }) {\n          throwBindingError('function ${humanName} called with ${
        arguments.length
      } arguments, expected ${argCount - 2} args!');\n        }`;
      if (needsDestructorStack) {
        invokerFnBody += "var destructors = [];\n";
      }
      var dtorStack = needsDestructorStack ? "destructors" : "null";
      var args1 = [
        "throwBindingError",
        "invoker",
        "fn",
        "runDestructors",
        "retType",
        "classParam",
      ];
      var args2 = [
        throwBindingError,
        cppInvokerFunc,
        cppTargetFunc,
        runDestructors,
        argTypes[0],
        argTypes[1],
      ];
      if (isClassMethodFunc) {
        invokerFnBody +=
          "var thisWired = classParam.toWireType(" + dtorStack + ", this);\n";
      }
      for (var i = 0; i < argCount - 2; ++i) {
        invokerFnBody +=
          "var arg" +
          i +
          "Wired = argType" +
          i +
          ".toWireType(" +
          dtorStack +
          ", arg" +
          i +
          "); // " +
          argTypes[i + 2].name +
          "\n";
        args1.push("argType" + i);
        args2.push(argTypes[i + 2]);
      }
      if (isClassMethodFunc) {
        argsListWired =
          "thisWired" + (argsListWired.length > 0 ? ", " : "") + argsListWired;
      }
      invokerFnBody +=
        (returns || isAsync ? "var rv = " : "") +
        "invoker(fn" +
        (argsListWired.length > 0 ? ", " : "") +
        argsListWired +
        ");\n";
      if (needsDestructorStack) {
        invokerFnBody += "runDestructors(destructors);\n";
      } else {
        for (var i = isClassMethodFunc ? 1 : 2; i < argTypes.length; ++i) {
          var paramName = i === 1 ? "thisWired" : "arg" + (i - 2) + "Wired";
          if (argTypes[i].destructorFunction !== null) {
            invokerFnBody +=
              paramName +
              "_dtor(" +
              paramName +
              "); // " +
              argTypes[i].name +
              "\n";
            args1.push(paramName + "_dtor");
            args2.push(argTypes[i].destructorFunction);
          }
        }
      }
      if (returns) {
        invokerFnBody +=
          "var ret = retType.fromWireType(rv);\n" + "return ret;\n";
      } else {
      }
      invokerFnBody += "}\n";
      args1.push(invokerFnBody);
      return newFunc(Function, args1).apply(null, args2);
    }
    function __embind_register_class_constructor(
      rawClassType,
      argCount,
      rawArgTypesAddr,
      invokerSignature,
      invoker,
      rawConstructor,
    ) {
      assert(argCount > 0);
      var rawArgTypes = heap32VectorToArray(argCount, rawArgTypesAddr);
      invoker = embind__requireFunction(invokerSignature, invoker);
      whenDependentTypesAreResolved([], [rawClassType], function (classType) {
        classType = classType[0];
        var humanName = `constructor ${classType.name}`;
        if (undefined === classType.registeredClass.constructor_body) {
          classType.registeredClass.constructor_body = [];
        }
        if (
          undefined !== classType.registeredClass.constructor_body[argCount - 1]
        ) {
          throw new BindingError(
            `Cannot register multiple constructors with identical number of parameters (${
              argCount - 1
            }) for class '${
              classType.name
            }'! Overload resolution is currently only performed using the parameter count, not actual type info!`,
          );
        }
        classType.registeredClass.constructor_body[argCount - 1] = () => {
          throwUnboundTypeError(
            `Cannot construct ${classType.name} due to unbound types`,
            rawArgTypes,
          );
        };
        whenDependentTypesAreResolved([], rawArgTypes, function (argTypes) {
          argTypes.splice(1, 0, null);
          classType.registeredClass.constructor_body[argCount - 1] =
            craftInvokerFunction(
              humanName,
              argTypes,
              null,
              invoker,
              rawConstructor,
            );
          return [];
        });
        return [];
      });
    }
    function __embind_register_class_function(
      rawClassType,
      methodName,
      argCount,
      rawArgTypesAddr,
      invokerSignature,
      rawInvoker,
      context,
      isPureVirtual,
      isAsync,
    ) {
      var rawArgTypes = heap32VectorToArray(argCount, rawArgTypesAddr);
      methodName = readLatin1String(methodName);
      rawInvoker = embind__requireFunction(invokerSignature, rawInvoker);
      whenDependentTypesAreResolved([], [rawClassType], function (classType) {
        classType = classType[0];
        var humanName = `${classType.name}.${methodName}`;
        if (methodName.startsWith("@@")) {
          methodName = Symbol[methodName.substring(2)];
        }
        if (isPureVirtual) {
          classType.registeredClass.pureVirtualFunctions.push(methodName);
        }
        function unboundTypesHandler() {
          throwUnboundTypeError(
            `Cannot call ${humanName} due to unbound types`,
            rawArgTypes,
          );
        }
        var proto = classType.registeredClass.instancePrototype;
        var method = proto[methodName];
        if (
          undefined === method ||
          (undefined === method.overloadTable &&
            method.className !== classType.name &&
            method.argCount === argCount - 2)
        ) {
          unboundTypesHandler.argCount = argCount - 2;
          unboundTypesHandler.className = classType.name;
          proto[methodName] = unboundTypesHandler;
        } else {
          ensureOverloadTable(proto, methodName, humanName);
          proto[methodName].overloadTable[argCount - 2] = unboundTypesHandler;
        }
        whenDependentTypesAreResolved([], rawArgTypes, function (argTypes) {
          var memberFunction = craftInvokerFunction(
            humanName,
            argTypes,
            classType,
            rawInvoker,
            context,
            isAsync,
          );
          if (undefined === proto[methodName].overloadTable) {
            memberFunction.argCount = argCount - 2;
            proto[methodName] = memberFunction;
          } else {
            proto[methodName].overloadTable[argCount - 2] = memberFunction;
          }
          return [];
        });
        return [];
      });
    }
    function handleAllocatorInit() {
      Object.assign(HandleAllocator.prototype, {
        get(id) {
          assert(this.allocated[id] !== undefined, `invalid handle: ${id}`);
          return this.allocated[id];
        },
        has(id) {
          return this.allocated[id] !== undefined;
        },
        allocate(handle) {
          var id = this.freelist.pop() || this.allocated.length;
          this.allocated[id] = handle;
          return id;
        },
        free(id) {
          assert(this.allocated[id] !== undefined);
          this.allocated[id] = undefined;
          this.freelist.push(id);
        },
      });
    }
    function HandleAllocator() {
      this.allocated = [undefined];
      this.freelist = [];
    }
    var emval_handles = new HandleAllocator();
    function __emval_decref(handle) {
      if (
        handle >= emval_handles.reserved &&
        0 === --emval_handles.get(handle).refcount
      ) {
        emval_handles.free(handle);
      }
    }
    function count_emval_handles() {
      var count = 0;
      for (
        var i = emval_handles.reserved;
        i < emval_handles.allocated.length;
        ++i
      ) {
        if (emval_handles.allocated[i] !== undefined) {
          ++count;
        }
      }
      return count;
    }
    function init_emval() {
      emval_handles.allocated.push(
        { value: undefined },
        { value: null },
        { value: true },
        { value: false },
      );
      emval_handles.reserved = emval_handles.allocated.length;
      Module["count_emval_handles"] = count_emval_handles;
    }
    var Emval = {
      toValue: (handle) => {
        if (!handle) {
          throwBindingError("Cannot use deleted val. handle = " + handle);
        }
        return emval_handles.get(handle).value;
      },
      toHandle: (value) => {
        switch (value) {
          case undefined:
            return 1;
          case null:
            return 2;
          case true:
            return 3;
          case false:
            return 4;
          default: {
            return emval_handles.allocate({ refcount: 1, value: value });
          }
        }
      },
    };
    function __embind_register_emval(rawType, name) {
      name = readLatin1String(name);
      registerType(rawType, {
        name: name,
        fromWireType: function (handle) {
          var rv = Emval.toValue(handle);
          __emval_decref(handle);
          return rv;
        },
        toWireType: function (destructors, value) {
          return Emval.toHandle(value);
        },
        argPackAdvance: 8,
        readValueFromPointer: simpleReadValueFromPointer,
        destructorFunction: null,
      });
    }
    function embindRepr(v) {
      if (v === null) {
        return "null";
      }
      var t = typeof v;
      if (t === "object" || t === "array" || t === "function") {
        return v.toString();
      } else {
        return "" + v;
      }
    }
    function floatReadValueFromPointer(name, shift) {
      switch (shift) {
        case 2:
          return function (pointer) {
            return this["fromWireType"](HEAPF32[pointer >> 2]);
          };
        case 3:
          return function (pointer) {
            return this["fromWireType"](HEAPF64[pointer >> 3]);
          };
        default:
          throw new TypeError("Unknown float type: " + name);
      }
    }
    function __embind_register_float(rawType, name, size) {
      var shift = getShiftFromSize(size);
      name = readLatin1String(name);
      registerType(rawType, {
        name: name,
        fromWireType: function (value) {
          return value;
        },
        toWireType: function (destructors, value) {
          if (typeof value != "number" && typeof value != "boolean") {
            throw new TypeError(
              `Cannot convert ${embindRepr(value)} to ${this.name}`,
            );
          }
          return value;
        },
        argPackAdvance: 8,
        readValueFromPointer: floatReadValueFromPointer(name, shift),
        destructorFunction: null,
      });
    }
    function __embind_register_function(
      name,
      argCount,
      rawArgTypesAddr,
      signature,
      rawInvoker,
      fn,
      isAsync,
    ) {
      var argTypes = heap32VectorToArray(argCount, rawArgTypesAddr);
      name = readLatin1String(name);
      rawInvoker = embind__requireFunction(signature, rawInvoker);
      exposePublicSymbol(
        name,
        function () {
          throwUnboundTypeError(
            `Cannot call ${name} due to unbound types`,
            argTypes,
          );
        },
        argCount - 1,
      );
      whenDependentTypesAreResolved([], argTypes, function (argTypes) {
        var invokerArgsArray = [argTypes[0], null].concat(argTypes.slice(1));
        replacePublicSymbol(
          name,
          craftInvokerFunction(
            name,
            invokerArgsArray,
            null,
            rawInvoker,
            fn,
            isAsync,
          ),
          argCount - 1,
        );
        return [];
      });
    }
    function integerReadValueFromPointer(name, shift, signed) {
      switch (shift) {
        case 0:
          return signed
            ? function readS8FromPointer(pointer) {
                return HEAP8[pointer];
              }
            : function readU8FromPointer(pointer) {
                return HEAPU8[pointer];
              };
        case 1:
          return signed
            ? function readS16FromPointer(pointer) {
                return HEAP16[pointer >> 1];
              }
            : function readU16FromPointer(pointer) {
                return HEAPU16[pointer >> 1];
              };
        case 2:
          return signed
            ? function readS32FromPointer(pointer) {
                return HEAP32[pointer >> 2];
              }
            : function readU32FromPointer(pointer) {
                return HEAPU32[pointer >> 2];
              };
        default:
          throw new TypeError("Unknown integer type: " + name);
      }
    }
    function __embind_register_integer(
      primitiveType,
      name,
      size,
      minRange,
      maxRange,
    ) {
      name = readLatin1String(name);
      if (maxRange === -1) {
        maxRange = 4294967295;
      }
      var shift = getShiftFromSize(size);
      var fromWireType = (value) => value;
      if (minRange === 0) {
        var bitshift = 32 - 8 * size;
        fromWireType = (value) => (value << bitshift) >>> bitshift;
      }
      var isUnsignedType = name.includes("unsigned");
      var checkAssertions = (value, toTypeName) => {
        if (typeof value != "number" && typeof value != "boolean") {
          throw new TypeError(
            `Cannot convert "${embindRepr(value)}" to ${toTypeName}`,
          );
        }
        if (value < minRange || value > maxRange) {
          throw new TypeError(
            `Passing a number "${embindRepr(
              value,
            )}" from JS side to C/C++ side to an argument of type "${name}", which is outside the valid range [${minRange}, ${maxRange}]!`,
          );
        }
      };
      var toWireType;
      if (isUnsignedType) {
        toWireType = function (destructors, value) {
          checkAssertions(value, this.name);
          return value >>> 0;
        };
      } else {
        toWireType = function (destructors, value) {
          checkAssertions(value, this.name);
          return value;
        };
      }
      registerType(primitiveType, {
        name: name,
        fromWireType: fromWireType,
        toWireType: toWireType,
        argPackAdvance: 8,
        readValueFromPointer: integerReadValueFromPointer(
          name,
          shift,
          minRange !== 0,
        ),
        destructorFunction: null,
      });
    }
    function __embind_register_memory_view(rawType, dataTypeIndex, name) {
      var typeMapping = [
        Int8Array,
        Uint8Array,
        Int16Array,
        Uint16Array,
        Int32Array,
        Uint32Array,
        Float32Array,
        Float64Array,
      ];
      var TA = typeMapping[dataTypeIndex];
      function decodeMemoryView(handle) {
        handle = handle >> 2;
        var heap = HEAPU32;
        var size = heap[handle];
        var data = heap[handle + 1];
        return new TA(heap.buffer, data, size);
      }
      name = readLatin1String(name);
      registerType(
        rawType,
        {
          name: name,
          fromWireType: decodeMemoryView,
          argPackAdvance: 8,
          readValueFromPointer: decodeMemoryView,
        },
        { ignoreDuplicateRegistrations: true },
      );
    }
    function __embind_register_std_string(rawType, name) {
      name = readLatin1String(name);
      var stdStringIsUTF8 = name === "std::string";
      registerType(rawType, {
        name: name,
        fromWireType: function (value) {
          var length = HEAPU32[value >> 2];
          var payload = value + 4;
          var str;
          if (stdStringIsUTF8) {
            var decodeStartPtr = payload;
            for (var i = 0; i <= length; ++i) {
              var currentBytePtr = payload + i;
              if (i == length || HEAPU8[currentBytePtr] == 0) {
                var maxRead = currentBytePtr - decodeStartPtr;
                var stringSegment = UTF8ToString(decodeStartPtr, maxRead);
                if (str === undefined) {
                  str = stringSegment;
                } else {
                  str += String.fromCharCode(0);
                  str += stringSegment;
                }
                decodeStartPtr = currentBytePtr + 1;
              }
            }
          } else {
            var a = new Array(length);
            for (var i = 0; i < length; ++i) {
              a[i] = String.fromCharCode(HEAPU8[payload + i]);
            }
            str = a.join("");
          }
          _free(value);
          return str;
        },
        toWireType: function (destructors, value) {
          if (value instanceof ArrayBuffer) {
            value = new Uint8Array(value);
          }
          var length;
          var valueIsOfTypeString = typeof value == "string";
          if (
            !(
              valueIsOfTypeString ||
              value instanceof Uint8Array ||
              value instanceof Uint8ClampedArray ||
              value instanceof Int8Array
            )
          ) {
            throwBindingError("Cannot pass non-string to std::string");
          }
          if (stdStringIsUTF8 && valueIsOfTypeString) {
            length = lengthBytesUTF8(value);
          } else {
            length = value.length;
          }
          var base = _malloc(4 + length + 1);
          var ptr = base + 4;
          HEAPU32[base >> 2] = length;
          checkInt32(length);
          if (stdStringIsUTF8 && valueIsOfTypeString) {
            stringToUTF8(value, ptr, length + 1);
          } else {
            if (valueIsOfTypeString) {
              for (var i = 0; i < length; ++i) {
                var charCode = value.charCodeAt(i);
                if (charCode > 255) {
                  _free(ptr);
                  throwBindingError(
                    "String has UTF-16 code units that do not fit in 8 bits",
                  );
                }
                HEAPU8[ptr + i] = charCode;
              }
            } else {
              for (var i = 0; i < length; ++i) {
                HEAPU8[ptr + i] = value[i];
              }
            }
          }
          if (destructors !== null) {
            destructors.push(_free, base);
          }
          return base;
        },
        argPackAdvance: 8,
        readValueFromPointer: simpleReadValueFromPointer,
        destructorFunction: function (ptr) {
          _free(ptr);
        },
      });
    }
    var UTF16Decoder =
      typeof TextDecoder != "undefined"
        ? new TextDecoder("utf-16le")
        : undefined;
    var UTF16ToString = (ptr, maxBytesToRead) => {
      assert(
        ptr % 2 == 0,
        "Pointer passed to UTF16ToString must be aligned to two bytes!",
      );
      var endPtr = ptr;
      var idx = endPtr >> 1;
      var maxIdx = idx + maxBytesToRead / 2;
      while (!(idx >= maxIdx) && HEAPU16[idx]) ++idx;
      endPtr = idx << 1;
      if (endPtr - ptr > 32 && UTF16Decoder)
        return UTF16Decoder.decode(HEAPU8.subarray(ptr, endPtr));
      var str = "";
      for (var i = 0; !(i >= maxBytesToRead / 2); ++i) {
        var codeUnit = HEAP16[(ptr + i * 2) >> 1];
        if (codeUnit == 0) break;
        str += String.fromCharCode(codeUnit);
      }
      return str;
    };
    var stringToUTF16 = (str, outPtr, maxBytesToWrite) => {
      assert(
        outPtr % 2 == 0,
        "Pointer passed to stringToUTF16 must be aligned to two bytes!",
      );
      assert(
        typeof maxBytesToWrite == "number",
        "stringToUTF16(str, outPtr, maxBytesToWrite) is missing the third parameter that specifies the length of the output buffer!",
      );
      if (maxBytesToWrite === undefined) {
        maxBytesToWrite = 2147483647;
      }
      if (maxBytesToWrite < 2) return 0;
      maxBytesToWrite -= 2;
      var startPtr = outPtr;
      var numCharsToWrite =
        maxBytesToWrite < str.length * 2 ? maxBytesToWrite / 2 : str.length;
      for (var i = 0; i < numCharsToWrite; ++i) {
        var codeUnit = str.charCodeAt(i);
        HEAP16[outPtr >> 1] = codeUnit;
        checkInt16(codeUnit);
        outPtr += 2;
      }
      HEAP16[outPtr >> 1] = 0;
      checkInt16(0);
      return outPtr - startPtr;
    };
    var lengthBytesUTF16 = (str) => str.length * 2;
    var UTF32ToString = (ptr, maxBytesToRead) => {
      assert(
        ptr % 4 == 0,
        "Pointer passed to UTF32ToString must be aligned to four bytes!",
      );
      var i = 0;
      var str = "";
      while (!(i >= maxBytesToRead / 4)) {
        var utf32 = HEAP32[(ptr + i * 4) >> 2];
        if (utf32 == 0) break;
        ++i;
        if (utf32 >= 65536) {
          var ch = utf32 - 65536;
          str += String.fromCharCode(55296 | (ch >> 10), 56320 | (ch & 1023));
        } else {
          str += String.fromCharCode(utf32);
        }
      }
      return str;
    };
    var stringToUTF32 = (str, outPtr, maxBytesToWrite) => {
      assert(
        outPtr % 4 == 0,
        "Pointer passed to stringToUTF32 must be aligned to four bytes!",
      );
      assert(
        typeof maxBytesToWrite == "number",
        "stringToUTF32(str, outPtr, maxBytesToWrite) is missing the third parameter that specifies the length of the output buffer!",
      );
      if (maxBytesToWrite === undefined) {
        maxBytesToWrite = 2147483647;
      }
      if (maxBytesToWrite < 4) return 0;
      var startPtr = outPtr;
      var endPtr = startPtr + maxBytesToWrite - 4;
      for (var i = 0; i < str.length; ++i) {
        var codeUnit = str.charCodeAt(i);
        if (codeUnit >= 55296 && codeUnit <= 57343) {
          var trailSurrogate = str.charCodeAt(++i);
          codeUnit =
            (65536 + ((codeUnit & 1023) << 10)) | (trailSurrogate & 1023);
        }
        HEAP32[outPtr >> 2] = codeUnit;
        checkInt32(codeUnit);
        outPtr += 4;
        if (outPtr + 4 > endPtr) break;
      }
      HEAP32[outPtr >> 2] = 0;
      checkInt32(0);
      return outPtr - startPtr;
    };
    var lengthBytesUTF32 = (str) => {
      var len = 0;
      for (var i = 0; i < str.length; ++i) {
        var codeUnit = str.charCodeAt(i);
        if (codeUnit >= 55296 && codeUnit <= 57343) ++i;
        len += 4;
      }
      return len;
    };
    var __embind_register_std_wstring = function (rawType, charSize, name) {
      name = readLatin1String(name);
      var decodeString, encodeString, getHeap, lengthBytesUTF, shift;
      if (charSize === 2) {
        decodeString = UTF16ToString;
        encodeString = stringToUTF16;
        lengthBytesUTF = lengthBytesUTF16;
        getHeap = () => HEAPU16;
        shift = 1;
      } else if (charSize === 4) {
        decodeString = UTF32ToString;
        encodeString = stringToUTF32;
        lengthBytesUTF = lengthBytesUTF32;
        getHeap = () => HEAPU32;
        shift = 2;
      }
      registerType(rawType, {
        name: name,
        fromWireType: function (value) {
          var length = HEAPU32[value >> 2];
          var HEAP = getHeap();
          var str;
          var decodeStartPtr = value + 4;
          for (var i = 0; i <= length; ++i) {
            var currentBytePtr = value + 4 + i * charSize;
            if (i == length || HEAP[currentBytePtr >> shift] == 0) {
              var maxReadBytes = currentBytePtr - decodeStartPtr;
              var stringSegment = decodeString(decodeStartPtr, maxReadBytes);
              if (str === undefined) {
                str = stringSegment;
              } else {
                str += String.fromCharCode(0);
                str += stringSegment;
              }
              decodeStartPtr = currentBytePtr + charSize;
            }
          }
          _free(value);
          return str;
        },
        toWireType: function (destructors, value) {
          if (!(typeof value == "string")) {
            throwBindingError(
              `Cannot pass non-string to C++ string type ${name}`,
            );
          }
          var length = lengthBytesUTF(value);
          var ptr = _malloc(4 + length + charSize);
          HEAPU32[ptr >> 2] = length >> shift;
          encodeString(value, ptr + 4, length + charSize);
          if (destructors !== null) {
            destructors.push(_free, ptr);
          }
          return ptr;
        },
        argPackAdvance: 8,
        readValueFromPointer: simpleReadValueFromPointer,
        destructorFunction: function (ptr) {
          _free(ptr);
        },
      });
    };
    function __embind_register_value_object(
      rawType,
      name,
      constructorSignature,
      rawConstructor,
      destructorSignature,
      rawDestructor,
    ) {
      structRegistrations[rawType] = {
        name: readLatin1String(name),
        rawConstructor: embind__requireFunction(
          constructorSignature,
          rawConstructor,
        ),
        rawDestructor: embind__requireFunction(
          destructorSignature,
          rawDestructor,
        ),
        fields: [],
      };
    }
    function __embind_register_value_object_field(
      structType,
      fieldName,
      getterReturnType,
      getterSignature,
      getter,
      getterContext,
      setterArgumentType,
      setterSignature,
      setter,
      setterContext,
    ) {
      structRegistrations[structType].fields.push({
        fieldName: readLatin1String(fieldName),
        getterReturnType: getterReturnType,
        getter: embind__requireFunction(getterSignature, getter),
        getterContext: getterContext,
        setterArgumentType: setterArgumentType,
        setter: embind__requireFunction(setterSignature, setter),
        setterContext: setterContext,
      });
    }
    function __embind_register_void(rawType, name) {
      name = readLatin1String(name);
      registerType(rawType, {
        isVoid: true,
        name: name,
        argPackAdvance: 0,
        fromWireType: function () {
          return undefined;
        },
        toWireType: function (destructors, o) {
          return undefined;
        },
      });
    }
    function __emval_incref(handle) {
      if (handle > 4) {
        emval_handles.get(handle).refcount += 1;
      }
    }
    function requireRegisteredType(rawType, humanName) {
      var impl = registeredTypes[rawType];
      if (undefined === impl) {
        throwBindingError(
          humanName + " has unknown type " + getTypeName(rawType),
        );
      }
      return impl;
    }
    function __emval_take_value(type, arg) {
      type = requireRegisteredType(type, "_emval_take_value");
      var v = type["readValueFromPointer"](arg);
      return Emval.toHandle(v);
    }
    var _abort = () => {
      abort("native code called abort()");
    };
    var _emscripten_memcpy_big = (dest, src, num) =>
      HEAPU8.copyWithin(dest, src, src + num);
    var getHeapMax = () => 2147483648;
    var _emscripten_get_now;
    _emscripten_get_now = () => performance.now();
    var growMemory = (size) => {
      var b = wasmMemory.buffer;
      var pages = (size - b.byteLength + 65535) >>> 16;
      try {
        wasmMemory.grow(pages);
        updateMemoryViews();
        return 1;
      } catch (e) {
        err(
          `growMemory: Attempted to grow heap from ${b.byteLength} bytes to ${size} bytes, but got error: ${e}`,
        );
      }
    };
    var _emscripten_resize_heap = (requestedSize) => {
      var oldSize = HEAPU8.length;
      requestedSize >>>= 0;
      assert(requestedSize > oldSize);
      var maxHeapSize = getHeapMax();
      if (requestedSize > maxHeapSize) {
        err(
          `Cannot enlarge memory, asked to go up to ${requestedSize} bytes, but the limit is ${maxHeapSize} bytes!`,
        );
        return false;
      }
      var alignUp = (x, multiple) =>
        x + ((multiple - (x % multiple)) % multiple);
      for (var cutDown = 1; cutDown <= 4; cutDown *= 2) {
        var overGrownHeapSize = oldSize * (1 + 0.2 / cutDown);
        overGrownHeapSize = Math.min(
          overGrownHeapSize,
          requestedSize + 100663296,
        );
        var newSize = Math.min(
          maxHeapSize,
          alignUp(Math.max(requestedSize, overGrownHeapSize), 65536),
        );
        var t0 = _emscripten_get_now();
        var replacement = growMemory(newSize);
        var t1 = _emscripten_get_now();
        out(
          `Heap resize call from ${oldSize} to ${newSize} took ${
            t1 - t0
          } msecs. Success: ${!!replacement}`,
        );
        if (replacement) {
          return true;
        }
      }
      err(
        `Failed to grow the heap from ${oldSize} bytes to ${newSize} bytes, not enough memory!`,
      );
      return false;
    };
    var ENV = {};
    var getExecutableName = () => thisProgram || "./this.program";
    var getEnvStrings = () => {
      if (!getEnvStrings.strings) {
        var lang =
          (
            (typeof navigator == "object" &&
              navigator.languages &&
              navigator.languages[0]) ||
            "C"
          ).replace("-", "_") + ".UTF-8";
        var env = {
          USER: "web_user",
          LOGNAME: "web_user",
          PATH: "/",
          PWD: "/",
          HOME: "/home/web_user",
          LANG: lang,
          _: getExecutableName(),
        };
        for (var x in ENV) {
          if (ENV[x] === undefined) delete env[x];
          else env[x] = ENV[x];
        }
        var strings = [];
        for (var x in env) {
          strings.push(`${x}=${env[x]}`);
        }
        getEnvStrings.strings = strings;
      }
      return getEnvStrings.strings;
    };
    var stringToAscii = (str, buffer) => {
      for (var i = 0; i < str.length; ++i) {
        assert(str.charCodeAt(i) === (str.charCodeAt(i) & 255));
        HEAP8[buffer++ >> 0] = str.charCodeAt(i);
        checkInt8(str.charCodeAt(i));
      }
      HEAP8[buffer >> 0] = 0;
      checkInt8(0);
    };
    var _environ_get = (__environ, environ_buf) => {
      var bufSize = 0;
      getEnvStrings().forEach(function (string, i) {
        var ptr = environ_buf + bufSize;
        HEAPU32[(__environ + i * 4) >> 2] = ptr;
        checkInt32(ptr);
        stringToAscii(string, ptr);
        bufSize += string.length + 1;
      });
      return 0;
    };
    var _environ_sizes_get = (penviron_count, penviron_buf_size) => {
      var strings = getEnvStrings();
      HEAPU32[penviron_count >> 2] = strings.length;
      checkInt32(strings.length);
      var bufSize = 0;
      strings.forEach(function (string) {
        bufSize += string.length + 1;
      });
      HEAPU32[penviron_buf_size >> 2] = bufSize;
      checkInt32(bufSize);
      return 0;
    };
    function _fd_close(fd) {
      try {
        var stream = SYSCALLS.getStreamFromFD(fd);
        FS.close(stream);
        return 0;
      } catch (e) {
        if (typeof FS == "undefined" || !(e.name === "ErrnoError")) throw e;
        return e.errno;
      }
    }
    var doReadv = (stream, iov, iovcnt, offset) => {
      var ret = 0;
      for (var i = 0; i < iovcnt; i++) {
        var ptr = HEAPU32[iov >> 2];
        var len = HEAPU32[(iov + 4) >> 2];
        iov += 8;
        var curr = FS.read(stream, HEAP8, ptr, len, offset);
        if (curr < 0) return -1;
        ret += curr;
        if (curr < len) break;
        if (typeof offset !== "undefined") {
          offset += curr;
        }
      }
      return ret;
    };
    function _fd_read(fd, iov, iovcnt, pnum) {
      try {
        var stream = SYSCALLS.getStreamFromFD(fd);
        var num = doReadv(stream, iov, iovcnt);
        HEAPU32[pnum >> 2] = num;
        checkInt32(num);
        return 0;
      } catch (e) {
        if (typeof FS == "undefined" || !(e.name === "ErrnoError")) throw e;
        return e.errno;
      }
    }
    function convertI32PairToI53Checked(lo, hi) {
      assert(lo == lo >>> 0 || lo == (lo | 0));
      assert(hi === (hi | 0));
      return (hi + 2097152) >>> 0 < 4194305 - !!lo
        ? (lo >>> 0) + hi * 4294967296
        : NaN;
    }
    function _fd_seek(fd, offset_low, offset_high, whence, newOffset) {
      var offset = convertI32PairToI53Checked(offset_low, offset_high);
      try {
        if (isNaN(offset)) return 61;
        var stream = SYSCALLS.getStreamFromFD(fd);
        FS.llseek(stream, offset, whence);
        (tempI64 = [
          stream.position >>> 0,
          ((tempDouble = stream.position),
          +Math.abs(tempDouble) >= 1
            ? tempDouble > 0
              ? +Math.floor(tempDouble / 4294967296) >>> 0
              : ~~+Math.ceil(
                  (tempDouble - +(~~tempDouble >>> 0)) / 4294967296,
                ) >>> 0
            : 0),
        ]),
          (HEAP32[newOffset >> 2] = tempI64[0]),
          (HEAP32[(newOffset + 4) >> 2] = tempI64[1]);
        checkInt64(stream.position);
        if (stream.getdents && offset === 0 && whence === 0)
          stream.getdents = null;
        return 0;
      } catch (e) {
        if (typeof FS == "undefined" || !(e.name === "ErrnoError")) throw e;
        return e.errno;
      }
    }
    var doWritev = (stream, iov, iovcnt, offset) => {
      var ret = 0;
      for (var i = 0; i < iovcnt; i++) {
        var ptr = HEAPU32[iov >> 2];
        var len = HEAPU32[(iov + 4) >> 2];
        iov += 8;
        var curr = FS.write(stream, HEAP8, ptr, len, offset);
        if (curr < 0) return -1;
        ret += curr;
        if (typeof offset !== "undefined") {
          offset += curr;
        }
      }
      return ret;
    };
    function _fd_write(fd, iov, iovcnt, pnum) {
      try {
        var stream = SYSCALLS.getStreamFromFD(fd);
        var num = doWritev(stream, iov, iovcnt);
        HEAPU32[pnum >> 2] = num;
        checkInt32(num);
        return 0;
      } catch (e) {
        if (typeof FS == "undefined" || !(e.name === "ErrnoError")) throw e;
        return e.errno;
      }
    }
    var isLeapYear = (year) =>
      year % 4 === 0 && (year % 100 !== 0 || year % 400 === 0);
    var arraySum = (array, index) => {
      var sum = 0;
      for (var i = 0; i <= index; sum += array[i++]) {}
      return sum;
    };
    var MONTH_DAYS_LEAP = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    var MONTH_DAYS_REGULAR = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    var addDays = (date, days) => {
      var newDate = new Date(date.getTime());
      while (days > 0) {
        var leap = isLeapYear(newDate.getFullYear());
        var currentMonth = newDate.getMonth();
        var daysInCurrentMonth = (leap ? MONTH_DAYS_LEAP : MONTH_DAYS_REGULAR)[
          currentMonth
        ];
        if (days > daysInCurrentMonth - newDate.getDate()) {
          days -= daysInCurrentMonth - newDate.getDate() + 1;
          newDate.setDate(1);
          if (currentMonth < 11) {
            newDate.setMonth(currentMonth + 1);
          } else {
            newDate.setMonth(0);
            newDate.setFullYear(newDate.getFullYear() + 1);
          }
        } else {
          newDate.setDate(newDate.getDate() + days);
          return newDate;
        }
      }
      return newDate;
    };
    var writeArrayToMemory = (array, buffer) => {
      assert(
        array.length >= 0,
        "writeArrayToMemory array must have a length (should be an array or typed array)",
      );
      HEAP8.set(array, buffer);
    };
    var _strftime = (s, maxsize, format, tm) => {
      var tm_zone = HEAP32[(tm + 40) >> 2];
      var date = {
        tm_sec: HEAP32[tm >> 2],
        tm_min: HEAP32[(tm + 4) >> 2],
        tm_hour: HEAP32[(tm + 8) >> 2],
        tm_mday: HEAP32[(tm + 12) >> 2],
        tm_mon: HEAP32[(tm + 16) >> 2],
        tm_year: HEAP32[(tm + 20) >> 2],
        tm_wday: HEAP32[(tm + 24) >> 2],
        tm_yday: HEAP32[(tm + 28) >> 2],
        tm_isdst: HEAP32[(tm + 32) >> 2],
        tm_gmtoff: HEAP32[(tm + 36) >> 2],
        tm_zone: tm_zone ? UTF8ToString(tm_zone) : "",
      };
      var pattern = UTF8ToString(format);
      var EXPANSION_RULES_1 = {
        "%c": "%a %b %d %H:%M:%S %Y",
        "%D": "%m/%d/%y",
        "%F": "%Y-%m-%d",
        "%h": "%b",
        "%r": "%I:%M:%S %p",
        "%R": "%H:%M",
        "%T": "%H:%M:%S",
        "%x": "%m/%d/%y",
        "%X": "%H:%M:%S",
        "%Ec": "%c",
        "%EC": "%C",
        "%Ex": "%m/%d/%y",
        "%EX": "%H:%M:%S",
        "%Ey": "%y",
        "%EY": "%Y",
        "%Od": "%d",
        "%Oe": "%e",
        "%OH": "%H",
        "%OI": "%I",
        "%Om": "%m",
        "%OM": "%M",
        "%OS": "%S",
        "%Ou": "%u",
        "%OU": "%U",
        "%OV": "%V",
        "%Ow": "%w",
        "%OW": "%W",
        "%Oy": "%y",
      };
      for (var rule in EXPANSION_RULES_1) {
        pattern = pattern.replace(
          new RegExp(rule, "g"),
          EXPANSION_RULES_1[rule],
        );
      }
      var WEEKDAYS = [
        "Sunday",
        "Monday",
        "Tuesday",
        "Wednesday",
        "Thursday",
        "Friday",
        "Saturday",
      ];
      var MONTHS = [
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
      ];
      function leadingSomething(value, digits, character) {
        var str = typeof value == "number" ? value.toString() : value || "";
        while (str.length < digits) {
          str = character[0] + str;
        }
        return str;
      }
      function leadingNulls(value, digits) {
        return leadingSomething(value, digits, "0");
      }
      function compareByDay(date1, date2) {
        function sgn(value) {
          return value < 0 ? -1 : value > 0 ? 1 : 0;
        }
        var compare;
        if ((compare = sgn(date1.getFullYear() - date2.getFullYear())) === 0) {
          if ((compare = sgn(date1.getMonth() - date2.getMonth())) === 0) {
            compare = sgn(date1.getDate() - date2.getDate());
          }
        }
        return compare;
      }
      function getFirstWeekStartDate(janFourth) {
        switch (janFourth.getDay()) {
          case 0:
            return new Date(janFourth.getFullYear() - 1, 11, 29);
          case 1:
            return janFourth;
          case 2:
            return new Date(janFourth.getFullYear(), 0, 3);
          case 3:
            return new Date(janFourth.getFullYear(), 0, 2);
          case 4:
            return new Date(janFourth.getFullYear(), 0, 1);
          case 5:
            return new Date(janFourth.getFullYear() - 1, 11, 31);
          case 6:
            return new Date(janFourth.getFullYear() - 1, 11, 30);
        }
      }
      function getWeekBasedYear(date) {
        var thisDate = addDays(
          new Date(date.tm_year + 1900, 0, 1),
          date.tm_yday,
        );
        var janFourthThisYear = new Date(thisDate.getFullYear(), 0, 4);
        var janFourthNextYear = new Date(thisDate.getFullYear() + 1, 0, 4);
        var firstWeekStartThisYear = getFirstWeekStartDate(janFourthThisYear);
        var firstWeekStartNextYear = getFirstWeekStartDate(janFourthNextYear);
        if (compareByDay(firstWeekStartThisYear, thisDate) <= 0) {
          if (compareByDay(firstWeekStartNextYear, thisDate) <= 0) {
            return thisDate.getFullYear() + 1;
          }
          return thisDate.getFullYear();
        }
        return thisDate.getFullYear() - 1;
      }
      var EXPANSION_RULES_2 = {
        "%a": (date) => WEEKDAYS[date.tm_wday].substring(0, 3),
        "%A": (date) => WEEKDAYS[date.tm_wday],
        "%b": (date) => MONTHS[date.tm_mon].substring(0, 3),
        "%B": (date) => MONTHS[date.tm_mon],
        "%C": (date) => {
          var year = date.tm_year + 1900;
          return leadingNulls((year / 100) | 0, 2);
        },
        "%d": (date) => leadingNulls(date.tm_mday, 2),
        "%e": (date) => leadingSomething(date.tm_mday, 2, " "),
        "%g": (date) => getWeekBasedYear(date).toString().substring(2),
        "%G": (date) => getWeekBasedYear(date),
        "%H": (date) => leadingNulls(date.tm_hour, 2),
        "%I": (date) => {
          var twelveHour = date.tm_hour;
          if (twelveHour == 0) twelveHour = 12;
          else if (twelveHour > 12) twelveHour -= 12;
          return leadingNulls(twelveHour, 2);
        },
        "%j": (date) =>
          leadingNulls(
            date.tm_mday +
              arraySum(
                isLeapYear(date.tm_year + 1900)
                  ? MONTH_DAYS_LEAP
                  : MONTH_DAYS_REGULAR,
                date.tm_mon - 1,
              ),
            3,
          ),
        "%m": (date) => leadingNulls(date.tm_mon + 1, 2),
        "%M": (date) => leadingNulls(date.tm_min, 2),
        "%n": () => "\n",
        "%p": (date) => {
          if (date.tm_hour >= 0 && date.tm_hour < 12) {
            return "AM";
          }
          return "PM";
        },
        "%S": (date) => leadingNulls(date.tm_sec, 2),
        "%t": () => "\t",
        "%u": (date) => date.tm_wday || 7,
        "%U": (date) => {
          var days = date.tm_yday + 7 - date.tm_wday;
          return leadingNulls(Math.floor(days / 7), 2);
        },
        "%V": (date) => {
          var val = Math.floor(
            (date.tm_yday + 7 - ((date.tm_wday + 6) % 7)) / 7,
          );
          if ((date.tm_wday + 371 - date.tm_yday - 2) % 7 <= 2) {
            val++;
          }
          if (!val) {
            val = 52;
            var dec31 = (date.tm_wday + 7 - date.tm_yday - 1) % 7;
            if (
              dec31 == 4 ||
              (dec31 == 5 && isLeapYear((date.tm_year % 400) - 1))
            ) {
              val++;
            }
          } else if (val == 53) {
            var jan1 = (date.tm_wday + 371 - date.tm_yday) % 7;
            if (jan1 != 4 && (jan1 != 3 || !isLeapYear(date.tm_year))) val = 1;
          }
          return leadingNulls(val, 2);
        },
        "%w": (date) => date.tm_wday,
        "%W": (date) => {
          var days = date.tm_yday + 7 - ((date.tm_wday + 6) % 7);
          return leadingNulls(Math.floor(days / 7), 2);
        },
        "%y": (date) => (date.tm_year + 1900).toString().substring(2),
        "%Y": (date) => date.tm_year + 1900,
        "%z": (date) => {
          var off = date.tm_gmtoff;
          var ahead = off >= 0;
          off = Math.abs(off) / 60;
          off = (off / 60) * 100 + (off % 60);
          return (ahead ? "+" : "-") + String("0000" + off).slice(-4);
        },
        "%Z": (date) => date.tm_zone,
        "%%": () => "%",
      };
      pattern = pattern.replace(/%%/g, "\0\0");
      for (var rule in EXPANSION_RULES_2) {
        if (pattern.includes(rule)) {
          pattern = pattern.replace(
            new RegExp(rule, "g"),
            EXPANSION_RULES_2[rule](date),
          );
        }
      }
      pattern = pattern.replace(/\0\0/g, "%");
      var bytes = intArrayFromString(pattern, false);
      if (bytes.length > maxsize) {
        return 0;
      }
      writeArrayToMemory(bytes, s);
      return bytes.length - 1;
    };
    var _strftime_l = (s, maxsize, format, tm, loc) =>
      _strftime(s, maxsize, format, tm);
    var FSNode = function (parent, name, mode, rdev) {
      if (!parent) {
        parent = this;
      }
      this.parent = parent;
      this.mount = parent.mount;
      this.mounted = null;
      this.id = FS.nextInode++;
      this.name = name;
      this.mode = mode;
      this.node_ops = {};
      this.stream_ops = {};
      this.rdev = rdev;
    };
    var readMode = 292 | 73;
    var writeMode = 146;
    Object.defineProperties(FSNode.prototype, {
      read: {
        get: function () {
          return (this.mode & readMode) === readMode;
        },
        set: function (val) {
          val ? (this.mode |= readMode) : (this.mode &= ~readMode);
        },
      },
      write: {
        get: function () {
          return (this.mode & writeMode) === writeMode;
        },
        set: function (val) {
          val ? (this.mode |= writeMode) : (this.mode &= ~writeMode);
        },
      },
      isFolder: {
        get: function () {
          return FS.isDir(this.mode);
        },
      },
      isDevice: {
        get: function () {
          return FS.isChrdev(this.mode);
        },
      },
    });
    FS.FSNode = FSNode;
    FS.createPreloadedFile = FS_createPreloadedFile;
    FS.staticInit();
    Module["FS_createPath"] = FS.createPath;
    Module["FS_createDataFile"] = FS.createDataFile;
    Module["FS_createPreloadedFile"] = FS.createPreloadedFile;
    Module["FS_unlink"] = FS.unlink;
    Module["FS_createLazyFile"] = FS.createLazyFile;
    Module["FS_createDevice"] = FS.createDevice;
    if (ENVIRONMENT_IS_NODE) {
      NODEFS.staticInit();
    }
    ERRNO_CODES = {
      EPERM: 63,
      ENOENT: 44,
      ESRCH: 71,
      EINTR: 27,
      EIO: 29,
      ENXIO: 60,
      E2BIG: 1,
      ENOEXEC: 45,
      EBADF: 8,
      ECHILD: 12,
      EAGAIN: 6,
      EWOULDBLOCK: 6,
      ENOMEM: 48,
      EACCES: 2,
      EFAULT: 21,
      ENOTBLK: 105,
      EBUSY: 10,
      EEXIST: 20,
      EXDEV: 75,
      ENODEV: 43,
      ENOTDIR: 54,
      EISDIR: 31,
      EINVAL: 28,
      ENFILE: 41,
      EMFILE: 33,
      ENOTTY: 59,
      ETXTBSY: 74,
      EFBIG: 22,
      ENOSPC: 51,
      ESPIPE: 70,
      EROFS: 69,
      EMLINK: 34,
      EPIPE: 64,
      EDOM: 18,
      ERANGE: 68,
      ENOMSG: 49,
      EIDRM: 24,
      ECHRNG: 106,
      EL2NSYNC: 156,
      EL3HLT: 107,
      EL3RST: 108,
      ELNRNG: 109,
      EUNATCH: 110,
      ENOCSI: 111,
      EL2HLT: 112,
      EDEADLK: 16,
      ENOLCK: 46,
      EBADE: 113,
      EBADR: 114,
      EXFULL: 115,
      ENOANO: 104,
      EBADRQC: 103,
      EBADSLT: 102,
      EDEADLOCK: 16,
      EBFONT: 101,
      ENOSTR: 100,
      ENODATA: 116,
      ETIME: 117,
      ENOSR: 118,
      ENONET: 119,
      ENOPKG: 120,
      EREMOTE: 121,
      ENOLINK: 47,
      EADV: 122,
      ESRMNT: 123,
      ECOMM: 124,
      EPROTO: 65,
      EMULTIHOP: 36,
      EDOTDOT: 125,
      EBADMSG: 9,
      ENOTUNIQ: 126,
      EBADFD: 127,
      EREMCHG: 128,
      ELIBACC: 129,
      ELIBBAD: 130,
      ELIBSCN: 131,
      ELIBMAX: 132,
      ELIBEXEC: 133,
      ENOSYS: 52,
      ENOTEMPTY: 55,
      ENAMETOOLONG: 37,
      ELOOP: 32,
      EOPNOTSUPP: 138,
      EPFNOSUPPORT: 139,
      ECONNRESET: 15,
      ENOBUFS: 42,
      EAFNOSUPPORT: 5,
      EPROTOTYPE: 67,
      ENOTSOCK: 57,
      ENOPROTOOPT: 50,
      ESHUTDOWN: 140,
      ECONNREFUSED: 14,
      EADDRINUSE: 3,
      ECONNABORTED: 13,
      ENETUNREACH: 40,
      ENETDOWN: 38,
      ETIMEDOUT: 73,
      EHOSTDOWN: 142,
      EHOSTUNREACH: 23,
      EINPROGRESS: 26,
      EALREADY: 7,
      EDESTADDRREQ: 17,
      EMSGSIZE: 35,
      EPROTONOSUPPORT: 66,
      ESOCKTNOSUPPORT: 137,
      EADDRNOTAVAIL: 4,
      ENETRESET: 39,
      EISCONN: 30,
      ENOTCONN: 53,
      ETOOMANYREFS: 141,
      EUSERS: 136,
      EDQUOT: 19,
      ESTALE: 72,
      ENOTSUP: 138,
      ENOMEDIUM: 148,
      EILSEQ: 25,
      EOVERFLOW: 61,
      ECANCELED: 11,
      ENOTRECOVERABLE: 56,
      EOWNERDEAD: 62,
      ESTRPIPE: 135,
    };
    InternalError = Module["InternalError"] = class InternalError extends (
      Error
    ) {
      constructor(message) {
        super(message);
        this.name = "InternalError";
      }
    };
    embind_init_charCodes();
    BindingError = Module["BindingError"] = class BindingError extends Error {
      constructor(message) {
        super(message);
        this.name = "BindingError";
      }
    };
    init_ClassHandle();
    init_embind();
    init_RegisteredPointer();
    UnboundTypeError = Module["UnboundTypeError"] = extendError(
      Error,
      "UnboundTypeError",
    );
    handleAllocatorInit();
    init_emval();
    function checkIncomingModuleAPI() {
      ignoredModuleProp("fetchSettings");
    }
    var wasmImports = {
      __assert_fail: ___assert_fail,
      __handle_stack_overflow: ___handle_stack_overflow,
      __syscall_fcntl64: ___syscall_fcntl64,
      __syscall_ioctl: ___syscall_ioctl,
      __syscall_openat: ___syscall_openat,
      __throw_exception_with_stack_trace: ___throw_exception_with_stack_trace,
      _embind_finalize_value_object: __embind_finalize_value_object,
      _embind_register_bigint: __embind_register_bigint,
      _embind_register_bool: __embind_register_bool,
      _embind_register_class: __embind_register_class,
      _embind_register_class_constructor: __embind_register_class_constructor,
      _embind_register_class_function: __embind_register_class_function,
      _embind_register_emval: __embind_register_emval,
      _embind_register_float: __embind_register_float,
      _embind_register_function: __embind_register_function,
      _embind_register_integer: __embind_register_integer,
      _embind_register_memory_view: __embind_register_memory_view,
      _embind_register_std_string: __embind_register_std_string,
      _embind_register_std_wstring: __embind_register_std_wstring,
      _embind_register_value_object: __embind_register_value_object,
      _embind_register_value_object_field: __embind_register_value_object_field,
      _embind_register_void: __embind_register_void,
      _emval_decref: __emval_decref,
      _emval_incref: __emval_incref,
      _emval_take_value: __emval_take_value,
      abort: _abort,
      emscripten_memcpy_big: _emscripten_memcpy_big,
      emscripten_resize_heap: _emscripten_resize_heap,
      environ_get: _environ_get,
      environ_sizes_get: _environ_sizes_get,
      fd_close: _fd_close,
      fd_read: _fd_read,
      fd_seek: _fd_seek,
      fd_write: _fd_write,
      strftime_l: _strftime_l,
    };
    var asm = createWasm();
    var ___wasm_call_ctors = createExportWrapper("__wasm_call_ctors");
    var _malloc = createExportWrapper("malloc");
    var _fflush = (Module["_fflush"] = createExportWrapper("fflush"));
    var _free = (Module["_free"] = createExportWrapper("free"));
    var ___errno_location = createExportWrapper("__errno_location");
    var ___getTypeName = createExportWrapper("__getTypeName");
    var __embind_initialize_bindings = (Module["__embind_initialize_bindings"] =
      createExportWrapper("_embind_initialize_bindings"));
    var _emscripten_builtin_memalign = createExportWrapper(
      "emscripten_builtin_memalign",
    );
    var ___trap = function () {
      return (___trap = Module["asm"]["__trap"]).apply(null, arguments);
    };
    var setTempRet0 = createExportWrapper("setTempRet0");
    var _emscripten_stack_init = function () {
      return (_emscripten_stack_init =
        Module["asm"]["emscripten_stack_init"]).apply(null, arguments);
    };
    var _emscripten_stack_get_free = function () {
      return (_emscripten_stack_get_free =
        Module["asm"]["emscripten_stack_get_free"]).apply(null, arguments);
    };
    var _emscripten_stack_get_base = function () {
      return (_emscripten_stack_get_base =
        Module["asm"]["emscripten_stack_get_base"]).apply(null, arguments);
    };
    var _emscripten_stack_get_end = function () {
      return (_emscripten_stack_get_end =
        Module["asm"]["emscripten_stack_get_end"]).apply(null, arguments);
    };
    var stackSave = createExportWrapper("stackSave");
    var stackRestore = createExportWrapper("stackRestore");
    var stackAlloc = createExportWrapper("stackAlloc");
    var _emscripten_stack_get_current = function () {
      return (_emscripten_stack_get_current =
        Module["asm"]["emscripten_stack_get_current"]).apply(null, arguments);
    };
    var ___cxa_decrement_exception_refcount = (Module[
      "___cxa_decrement_exception_refcount"
    ] = createExportWrapper("__cxa_decrement_exception_refcount"));
    var ___cxa_increment_exception_refcount = (Module[
      "___cxa_increment_exception_refcount"
    ] = createExportWrapper("__cxa_increment_exception_refcount"));
    var ___cxa_demangle = createExportWrapper("__cxa_demangle");
    var ___thrown_object_from_unwind_exception = (Module[
      "___thrown_object_from_unwind_exception"
    ] = createExportWrapper("__thrown_object_from_unwind_exception"));
    var ___get_exception_message = (Module["___get_exception_message"] =
      createExportWrapper("__get_exception_message"));
    var ___set_stack_limits = (Module["___set_stack_limits"] =
      createExportWrapper("__set_stack_limits"));
    var dynCall_jiji = (Module["dynCall_jiji"] =
      createExportWrapper("dynCall_jiji"));
    var dynCall_viijii = (Module["dynCall_viijii"] =
      createExportWrapper("dynCall_viijii"));
    var dynCall_iiiiij = (Module["dynCall_iiiiij"] =
      createExportWrapper("dynCall_iiiiij"));
    var dynCall_iiiiijj = (Module["dynCall_iiiiijj"] =
      createExportWrapper("dynCall_iiiiijj"));
    var dynCall_iiiiiijj = (Module["dynCall_iiiiiijj"] =
      createExportWrapper("dynCall_iiiiiijj"));
    function intArrayFromBase64(s) {
      if (typeof ENVIRONMENT_IS_NODE != "undefined" && ENVIRONMENT_IS_NODE) {
        var buf = Buffer.from(s, "base64");
        return new Uint8Array(
          buf["buffer"],
          buf["byteOffset"],
          buf["byteLength"],
        );
      }
      try {
        var decoded = atob(s);
        var bytes = new Uint8Array(decoded.length);
        for (var i = 0; i < decoded.length; ++i) {
          bytes[i] = decoded.charCodeAt(i);
        }
        return bytes;
      } catch (_) {
        throw new Error("Converting base64 string to bytes failed.");
      }
    }
    Module["addRunDependency"] = addRunDependency;
    Module["removeRunDependency"] = removeRunDependency;
    Module["FS_createPath"] = FS.createPath;
    Module["FS_createDataFile"] = FS.createDataFile;
    Module["FS_createLazyFile"] = FS.createLazyFile;
    Module["FS_createDevice"] = FS.createDevice;
    Module["FS_unlink"] = FS.unlink;
    Module["FS_createPreloadedFile"] = FS.createPreloadedFile;
    Module["FS"] = FS;
    var missingLibrarySymbols = [
      "writeI53ToI64",
      "writeI53ToI64Clamped",
      "writeI53ToI64Signaling",
      "writeI53ToU64Clamped",
      "writeI53ToU64Signaling",
      "readI53FromI64",
      "readI53FromU64",
      "convertI32PairToI53",
      "convertU32PairToI53",
      "exitJS",
      "ydayFromDate",
      "inetPton4",
      "inetNtop4",
      "inetPton6",
      "inetNtop6",
      "readSockaddr",
      "writeSockaddr",
      "getHostByName",
      "traverseStack",
      "getCallstack",
      "emscriptenLog",
      "convertPCtoSourceLocation",
      "readEmAsmArgs",
      "jstoi_q",
      "jstoi_s",
      "listenOnce",
      "autoResumeAudioContext",
      "handleException",
      "runtimeKeepalivePush",
      "runtimeKeepalivePop",
      "callUserCallback",
      "maybeExit",
      "safeSetTimeout",
      "asmjsMangle",
      "getNativeTypeSize",
      "STACK_SIZE",
      "STACK_ALIGN",
      "POINTER_SIZE",
      "ASSERTIONS",
      "getCFunc",
      "ccall",
      "cwrap",
      "uleb128Encode",
      "sigToWasmTypes",
      "generateFuncType",
      "convertJsFunctionToWasm",
      "getEmptyTableSlot",
      "updateTableMap",
      "getFunctionAddress",
      "addFunction",
      "removeFunction",
      "reallyNegative",
      "unSign",
      "strLen",
      "reSign",
      "formatString",
      "intArrayToString",
      "AsciiToString",
      "stringToNewUTF8",
      "registerKeyEventCallback",
      "maybeCStringToJsString",
      "findEventTarget",
      "findCanvasEventTarget",
      "getBoundingClientRect",
      "fillMouseEventData",
      "registerMouseEventCallback",
      "registerWheelEventCallback",
      "registerUiEventCallback",
      "registerFocusEventCallback",
      "fillDeviceOrientationEventData",
      "registerDeviceOrientationEventCallback",
      "fillDeviceMotionEventData",
      "registerDeviceMotionEventCallback",
      "screenOrientation",
      "fillOrientationChangeEventData",
      "registerOrientationChangeEventCallback",
      "fillFullscreenChangeEventData",
      "registerFullscreenChangeEventCallback",
      "JSEvents_requestFullscreen",
      "JSEvents_resizeCanvasForFullscreen",
      "registerRestoreOldStyle",
      "hideEverythingExceptGivenElement",
      "restoreHiddenElements",
      "setLetterbox",
      "softFullscreenResizeWebGLRenderTarget",
      "doRequestFullscreen",
      "fillPointerlockChangeEventData",
      "registerPointerlockChangeEventCallback",
      "registerPointerlockErrorEventCallback",
      "requestPointerLock",
      "fillVisibilityChangeEventData",
      "registerVisibilityChangeEventCallback",
      "registerTouchEventCallback",
      "fillGamepadEventData",
      "registerGamepadEventCallback",
      "registerBeforeUnloadEventCallback",
      "fillBatteryEventData",
      "battery",
      "registerBatteryEventCallback",
      "setCanvasElementSize",
      "getCanvasElementSize",
      "checkWasiClock",
      "wasiRightsToMuslOFlags",
      "wasiOFlagsToMuslOFlags",
      "createDyncallWrapper",
      "setImmediateWrapped",
      "clearImmediateWrapped",
      "polyfillSetImmediate",
      "getPromise",
      "makePromise",
      "idsToPromises",
      "makePromiseCallback",
      "setMainLoop",
      "getSocketFromFD",
      "getSocketAddress",
      "_setNetworkCallback",
      "heapObjectForWebGLType",
      "heapAccessShiftForWebGLHeap",
      "webgl_enable_ANGLE_instanced_arrays",
      "webgl_enable_OES_vertex_array_object",
      "webgl_enable_WEBGL_draw_buffers",
      "webgl_enable_WEBGL_multi_draw",
      "emscriptenWebGLGet",
      "computeUnpackAlignedImageSize",
      "colorChannelsInGlTextureFormat",
      "emscriptenWebGLGetTexPixelData",
      "__glGenObject",
      "emscriptenWebGLGetUniform",
      "webglGetUniformLocation",
      "webglPrepareUniformLocationsBeforeFirstUse",
      "webglGetLeftBracePos",
      "emscriptenWebGLGetVertexAttrib",
      "__glGetActiveAttribOrUniform",
      "writeGLArray",
      "registerWebGlEventCallback",
      "runAndAbortIfError",
      "SDL_unicode",
      "SDL_ttfContext",
      "SDL_audio",
      "GLFW_Window",
      "ALLOC_NORMAL",
      "ALLOC_STACK",
      "allocate",
      "writeStringToMemory",
      "writeAsciiToMemory",
      "registerInheritedInstance",
      "unregisterInheritedInstance",
      "enumReadValueFromPointer",
      "validateThis",
      "getStringOrSymbol",
      "craftEmvalAllocator",
      "emval_get_global",
      "emval_lookupTypes",
      "emval_allocateDestructors",
      "emval_addMethodCaller",
    ];
    missingLibrarySymbols.forEach(missingLibrarySymbol);
    var unexportedSymbols = [
      "run",
      "addOnPreRun",
      "addOnInit",
      "addOnPreMain",
      "addOnExit",
      "addOnPostRun",
      "FS_createFolder",
      "FS_createLink",
      "out",
      "err",
      "callMain",
      "abort",
      "keepRuntimeAlive",
      "wasmMemory",
      "stackAlloc",
      "stackSave",
      "stackRestore",
      "getTempRet0",
      "setTempRet0",
      "writeStackCookie",
      "checkStackCookie",
      "convertI32PairToI53Checked",
      "ptrToString",
      "zeroMemory",
      "getHeapMax",
      "growMemory",
      "ENV",
      "setStackLimits",
      "MONTH_DAYS_REGULAR",
      "MONTH_DAYS_LEAP",
      "MONTH_DAYS_REGULAR_CUMULATIVE",
      "MONTH_DAYS_LEAP_CUMULATIVE",
      "isLeapYear",
      "arraySum",
      "addDays",
      "ERRNO_CODES",
      "ERRNO_MESSAGES",
      "setErrNo",
      "DNS",
      "Protocols",
      "Sockets",
      "initRandomFill",
      "randomFill",
      "timers",
      "warnOnce",
      "UNWIND_CACHE",
      "readEmAsmArgsArray",
      "getExecutableName",
      "dynCallLegacy",
      "getDynCaller",
      "dynCall",
      "asyncLoad",
      "alignMemory",
      "mmapAlloc",
      "handleAllocatorInit",
      "HandleAllocator",
      "freeTableIndexes",
      "functionsInTableMap",
      "setValue",
      "getValue",
      "PATH",
      "PATH_FS",
      "UTF8Decoder",
      "UTF8ArrayToString",
      "UTF8ToString",
      "stringToUTF8Array",
      "stringToUTF8",
      "lengthBytesUTF8",
      "intArrayFromString",
      "stringToAscii",
      "UTF16Decoder",
      "UTF16ToString",
      "stringToUTF16",
      "lengthBytesUTF16",
      "UTF32ToString",
      "stringToUTF32",
      "lengthBytesUTF32",
      "stringToUTF8OnStack",
      "writeArrayToMemory",
      "JSEvents",
      "specialHTMLTargets",
      "currentFullscreenStrategy",
      "restoreOldWindowedStyle",
      "demangle",
      "demangleAll",
      "jsStackTrace",
      "stackTrace",
      "ExitStatus",
      "getEnvStrings",
      "doReadv",
      "doWritev",
      "promiseMap",
      "getExceptionMessageCommon",
      "getCppExceptionTag",
      "getCppExceptionThrownObjectFromWebAssemblyException",
      "incrementExceptionRefcount",
      "decrementExceptionRefcount",
      "getExceptionMessage",
      "Browser",
      "wget",
      "SYSCALLS",
      "preloadPlugins",
      "FS_modeStringToFlags",
      "FS_getMode",
      "FS_stdin_getChar_buffer",
      "FS_stdin_getChar",
      "MEMFS",
      "TTY",
      "PIPEFS",
      "SOCKFS",
      "tempFixedLengthArray",
      "miniTempWebGLFloatBuffers",
      "miniTempWebGLIntBuffers",
      "GL",
      "emscripten_webgl_power_preferences",
      "AL",
      "GLUT",
      "EGL",
      "GLEW",
      "IDBStore",
      "SDL",
      "SDL_gfx",
      "GLFW",
      "allocateUTF8",
      "allocateUTF8OnStack",
      "InternalError",
      "BindingError",
      "throwInternalError",
      "throwBindingError",
      "registeredTypes",
      "awaitingDependencies",
      "typeDependencies",
      "tupleRegistrations",
      "structRegistrations",
      "sharedRegisterType",
      "whenDependentTypesAreResolved",
      "embind_charCodes",
      "embind_init_charCodes",
      "readLatin1String",
      "getTypeName",
      "heap32VectorToArray",
      "requireRegisteredType",
      "UnboundTypeError",
      "PureVirtualError",
      "init_embind",
      "throwUnboundTypeError",
      "ensureOverloadTable",
      "exposePublicSymbol",
      "replacePublicSymbol",
      "extendError",
      "createNamedFunction",
      "embindRepr",
      "registeredInstances",
      "getBasestPointer",
      "getInheritedInstance",
      "getInheritedInstanceCount",
      "getLiveInheritedInstances",
      "registeredPointers",
      "registerType",
      "getShiftFromSize",
      "integerReadValueFromPointer",
      "floatReadValueFromPointer",
      "simpleReadValueFromPointer",
      "runDestructors",
      "newFunc",
      "craftInvokerFunction",
      "embind__requireFunction",
      "genericPointerToWireType",
      "constNoSmartPtrRawPointerToWireType",
      "nonConstNoSmartPtrRawPointerToWireType",
      "init_RegisteredPointer",
      "RegisteredPointer",
      "RegisteredPointer_getPointee",
      "RegisteredPointer_destructor",
      "RegisteredPointer_deleteObject",
      "RegisteredPointer_fromWireType",
      "runDestructor",
      "releaseClassHandle",
      "finalizationRegistry",
      "detachFinalizer_deps",
      "detachFinalizer",
      "attachFinalizer",
      "makeClassHandle",
      "init_ClassHandle",
      "ClassHandle",
      "ClassHandle_isAliasOf",
      "throwInstanceAlreadyDeleted",
      "ClassHandle_clone",
      "ClassHandle_delete",
      "deletionQueue",
      "ClassHandle_isDeleted",
      "ClassHandle_deleteLater",
      "flushPendingDeletes",
      "delayFunction",
      "setDelayFunction",
      "RegisteredClass",
      "shallowCopyInternalPointer",
      "downcastPointer",
      "upcastPointer",
      "char_0",
      "char_9",
      "makeLegalFunctionName",
      "emval_handles",
      "emval_symbols",
      "init_emval",
      "count_emval_handles",
      "Emval",
      "emval_newers",
      "emval_methodCallers",
      "emval_registeredMethods",
      "NODEFS",
    ];
    unexportedSymbols.forEach(unexportedRuntimeSymbol);
    var calledRun;
    dependenciesFulfilled = function runCaller() {
      if (!calledRun) run();
      if (!calledRun) dependenciesFulfilled = runCaller;
    };
    function stackCheckInit() {
      _emscripten_stack_init();
      writeStackCookie();
    }
    function run() {
      if (runDependencies > 0) {
        return;
      }
      stackCheckInit();
      preRun();
      if (runDependencies > 0) {
        return;
      }
      function doRun() {
        if (calledRun) return;
        calledRun = true;
        Module["calledRun"] = true;
        if (ABORT) return;
        initRuntime();
        readyPromiseResolve(Module);
        if (Module["onRuntimeInitialized"]) Module["onRuntimeInitialized"]();
        assert(
          !Module["_main"],
          'compiled without a main, but one is present. if you added it from JS, use Module["onRuntimeInitialized"]',
        );
        postRun();
      }
      if (Module["setStatus"]) {
        Module["setStatus"]("Running...");
        setTimeout(function () {
          setTimeout(function () {
            Module["setStatus"]("");
          }, 1);
          doRun();
        }, 1);
      } else {
        doRun();
      }
      checkStackCookie();
    }
    if (Module["preInit"]) {
      if (typeof Module["preInit"] == "function")
        Module["preInit"] = [Module["preInit"]];
      while (Module["preInit"].length > 0) {
        Module["preInit"].pop()();
      }
    }
    run();

    return moduleArg.ready;
  };
})();
export default privateer_module;

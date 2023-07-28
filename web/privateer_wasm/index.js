console.log("Script running")
    privateer_module().then((module_) => {
        console.log("WASM LOADED");
                
        var button = document.getElementById("send-file-button");
        var fileInput = document.getElementById("upload");

        button.onclick = function () {
            var files = fileInput.files;
            var reader = new FileReader();
            reader.onload = function () {
                x = module_.read_structure(reader.result, files[0].name)
                let parent = document.getElementById("svg");

                for (var i = 0; i < x.size(); i++) {
                  const span = document.createElement('div')
                  span.innerHTML = x.get(i);
                  parent.appendChild(span);
                }


        };
        if(files[0]) {
            reader.readAsText(files[0]);
        }
        };
    
    
    })
const { contextBridge, ipcRenderer } = require("electron");

// Expose API to frontend — also provide legacy "vsnp" alias for compatibility
contextBridge.exposeInMainWorld("kraken", {
  selectPath: (opts) => ipcRenderer.invoke("select-path", opts),
  apiBase: process.env.VITE_API_URL || ""
});

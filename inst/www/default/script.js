/* jshint esversion:6 */

$(document).ready(() => {
    const observer = new MutationObserver(
        (mutationsList, observer) => {
            mutationsList.forEach((record) => {
                Array.from(record.addedNodes)
                    .filter(node => $(node).hasClass("array"))
                    .forEach((node) => {
                        const resolution = $(node).data("resolution");
                        const container = $("<div class='array-container'>")
                              .appendTo($(node).parent())
                        ;
                        const body = $("<div class='array-container-body'>")
                              .append($(node).detach())
                              .prependTo(container)
                        ;
                        const headingText = $(`<span>Resolution ${resolution}</span>`);
                        const headingIcon = $("<i>");
                        const heading = $("<div class='array-container-heading'>")
                              .click(() => {
                                  body.toggleClass("collapsed");
                                  if (body.hasClass("collapsed")) {
                                      headingIcon.attr("class", "glyphicon glyphicon-menu-right");
                                  } else {
                                      headingIcon.attr("class", "glyphicon glyphicon-menu-left");
                                  }
                              })
                              .append(headingIcon)
                              .append(headingText)
                              .prependTo(container)
                        ;
                        heading.click().click();
                    });
            });
        }
    );

    observer.observe(
        document.getElementById("array"),
        { childList: true }
    );
});

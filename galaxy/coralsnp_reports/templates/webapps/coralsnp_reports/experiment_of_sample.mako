<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Experiment of sample with affy id ${affy_id}
        </h3>
        <table align="center" class="colored">
            %if len(experiments) == 0:
                <tr>
                    <td colspan="2">
                        This sample has no experiment
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Seq Facility</td>
                    <td>Array Version</td>
                    <td>Result Folder Name</td>
                    <td>Plate Barcode</td>
                </tr>
                <% ctr = 0 %>
                %for experiment in experiments:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${experiment[0]}</td>
                        <td>${experiment[1]}</td>
                        <td>${experiment[2]}</td>
                        <td>${experiment[3]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

